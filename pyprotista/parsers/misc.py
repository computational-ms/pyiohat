"""Parser handler."""

import IsoSpecPy as iso
import regex as re
from IsoSpecPy.PeriodicTbl import symbol_to_masses


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
