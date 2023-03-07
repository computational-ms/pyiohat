"""Parser handler."""

import IsoSpecPy as iso
import numpy as np
import pandas as pd
import regex as re
from IsoSpecPy.PeriodicTbl import symbol_to_masses


def get_atom_counts(sequences, modifications, compositions):
    """Get the count occurrence of each atom in a given sequence with respective mods.

    Args:
        sequences (np.array): numpy string array with sequences
        modifications (np.array): numpy string array with modifications
        compositions (dict of dict): dict with each amino acid and modification composition

    Returns:
        elements (list): ordered uniquely occurring elements
        atom_counts (np.array): numpy int array with PSM-level atomic composition
    """
    unique_elements = set()
    for composition in compositions.values():
        for element in composition.keys():
            unique_elements.add(element)

    elements = ["C", "H"] + sorted(unique_elements - set(["C", "H"]))

    atom_counts = np.zeros(shape=(len(sequences), len(elements)), dtype=int)
    for aa_or_mod in compositions.keys():
        ordered_element_multiplier = []
        for element in elements:
            ordered_element_multiplier += [compositions[aa_or_mod].get(element, 0)]
        ordered_element_multiplier = np.array(ordered_element_multiplier)
        # Only amino acids are length 1
        if len(aa_or_mod) == 1:
            atom_counts += np.outer(
                np.char.count(sequences, aa_or_mod), ordered_element_multiplier
            )
        else:
            mod_counts = []
            escaped_mod_name = re.escape(aa_or_mod)
            search_pattern = re.compile(
                rf"(^{escaped_mod_name}:\d+)(?=;)|(?<=;)({escaped_mod_name}:\d+)(?=;)|(?<=;)({escaped_mod_name}:\d+$)|^({escaped_mod_name}:\d+)$"
            )
            for mod in modifications:
                mod_counts.append(
                    len(
                        re.findall(
                            search_pattern,
                            mod,
                        )
                    )
                )
            atom_counts += np.outer(mod_counts, ordered_element_multiplier)
    # Remove water (peptide bonds)
    water = np.zeros(shape=(1, len(elements)), dtype=int)
    water[0, elements.index("H")] = 2
    water[0, elements.index("O")] = 1
    atom_counts -= np.outer(np.maximum(np.char.str_len(sequences) - 1, 0), water)

    return elements, atom_counts


def get_compositions_and_monoisotopic_masses(
    sequences, modifications, compositions, isotopic_distributions
):
    """Get compositions and monoisotopic masses per PSM.

    Args:
        sequences (np.array): numpy string array with sequences
        modifications (np.array): numpy string array with modifications
        compositions (dict of dict): dict with each amino acid and modification composition
        isotopic_distributions (dict of dict): isotopic distributions as provided by chemical_composition

    Returns:
        chemical_compositions (pd.Series): PSM-level chemical compositions in hill notation
        monoisotopic_masses (np.array): numpy float array with PSM-level monoisotopic masses
    """
    # Get atom counts
    elements, atom_counts = get_atom_counts(
        sequences=sequences, modifications=modifications, compositions=compositions
    )
    for i, element in enumerate(elements):
        if i == 0:
            chemical_compositions = np.char.add(
                elements[i] + "(", atom_counts.astype(str)[:, i]
            )
        else:
            next_element_string = np.char.add(
                ")" + elements[i] + "(", atom_counts.astype(str)[:, i]
            )
            chemical_compositions = np.char.add(
                chemical_compositions, next_element_string
            )
        if i == len(elements) - 1:
            chemical_compositions = np.char.add(chemical_compositions, ")")
    # Remove empties
    chemical_compositions = pd.Series(chemical_compositions).str.replace(
        r"\d*[A-Z][a-z]?\(0\)", "", regex=True
    )

    isotope_mass_lookup = {}
    for element, isotope_data in isotopic_distributions.items():
        if isinstance(isotope_data, dict):
            # this has extra data
            isotope_list = []
            for fraction in isotope_data:
                isotope_list.extend(isotope_data[fraction])
        else:
            isotope_list = isotope_data

        for iso_abund in isotope_list:
            isotope_mass_key = f"{round(iso_abund[0])}{element}"
            isotope_mass_lookup[isotope_mass_key] = iso_abund[0]

    monoisotopic_masses = []
    for element in elements:
        try:
            mass = max(isotopic_distributions[element], key=lambda d: d[1])[0]
        except KeyError:
            # Element is an isotope
            mass = isotope_mass_lookup[element]
        monoisotopic_masses.extend([mass])
    monoisotopic_masses = np.array(monoisotopic_masses)
    monoisotopic_masses = (atom_counts * monoisotopic_masses).sum(axis=1)

    return chemical_compositions, monoisotopic_masses


def sort_mods(data):
    """Sort mods based on position in peptide and ensures unique mod:pos pairs.

    Args:
        data (list): list of modifications with style "Mod:position"
    Returns:
        sorted_formatted_mods (str): String with sorted mods in style "Mod1:pos1;Modn:posn"
    """

    def _sort_with_positional_priority(mod):
        if mod == "":
            return None
        else:
            name, _, pos = mod.rpartition(":")
            return int(pos), name

    sorted_formatted_mods = ";".join(
        sorted(set(data), key=_sort_with_positional_priority)
    )
    return sorted_formatted_mods


def init_custom_cc(function, proton):
    """Initialize function for multiprocessing by providing 'global' attribute.

    Args:
        function (function): function to be appended with cc attribute
        proton (float): proton mass
    """

    def _calc_mz(mass, charge):
        return (mass + (charge * proton)) / charge

    function.calc_mz = _calc_mz


def get_isotopologue_accuracy(composition, charge, exp_mz):
    """Compute the isotopologue accuracy.

    Only the accuracy of the isotopologue closest to the experimental mass is reported.

    Args:
        composition (str): chemical composition in hill notation
        charge (int): charge
        exp_mz (float): experimental spectrum mz
    Returns:
        isotopologue_acc (float): accuracy in ppm
    """
    atom_counts = None
    isotope_masses = None
    isotope_probs = None
    replaced_composition = composition
    static_isotopes = re.findall(r"(?<=\))(\d+)(\w+)(?:\()(\d+)", composition)
    if len(static_isotopes) != 0:
        atom_counts = []
        isotope_masses = []
        for isotope in static_isotopes:
            mass, element, number = isotope
            replaced_composition = replaced_composition.replace(
                f"{mass}{element}({number})", ""
            )
            isotope_masses.append(
                [
                    [
                        m
                        for m in symbol_to_masses[element]
                        if str(round(m)).startswith(mass)
                    ][0]
                ]
            )
            atom_counts.append(int(number))
        isotope_probs = len(atom_counts) * [[1.0]]
    if replaced_composition is not None:
        formula = replaced_composition
    else:
        formula = composition
    formula = formula.replace("(", "").replace(")", "")
    isotopologue_masses = list(
        iso.IsoThreshold(
            formula=formula,
            threshold=0.02,
            charge=1,
            get_confs=True,
            atomCounts=atom_counts,
            isotopeMasses=isotope_masses,
            isotopeProbabilities=isotope_probs,
        ).masses
    )
    isotopologue_mzs = [
        get_isotopologue_accuracy.calc_mz(mass, charge) for mass in isotopologue_masses
    ]
    # Report only most accurate mass
    isotopologue_mz = min(
        isotopologue_mzs,
        key=lambda x: abs(exp_mz - x),
    )
    isotopologue_acc = (exp_mz - isotopologue_mz) / isotopologue_mz * 1e6
    return isotopologue_acc
