"""Engine parser."""
import numpy as np
import pandas as pd
import regex as re
import xml.etree.ElementTree as etree
import warnings
from pyprotista.parsers.ident_base_parser import IdentBaseParser


def get_xml_data(xml_file, mapping_dict):
    """For loop over one xml file using xml.etree.ElementTree.iterparse.

    Provide temporary variables for iteration.
    Call functions get_peptide_lookup in stage 1, get_modification_mass_map in stage 2 and get_spec_records in stage 3.

    Args:
        xml_file (mzid): input file
        mapping_dict (dict)

    Returns:
        version (str): contains the version of the mzid
        peptide_lookup (dict): contains peptide id's with their sequences and modifications
        modification_mass_map (dict): contains name and massDelta of Modifications
        spec_records (list): list of dicts containing single_specs
        fixed_mods (dict): contains Modifications from modification_mass_map where fixedMod = true
    """
    version = ""
    peptide_lookup = {}
    modification_mass_map = {}
    spec_records = []
    fixed_mods = {}

    # Temporary variables, get overwritten multiple times during iteration
    stage = 0
    modifications = []
    sequence = {}
    spec_results = {}
    spec_ident_items = []
    mod_name = ""

    for event, entry in etree.iterparse(xml_file):
        entry_tag = entry.tag
        if entry_tag.endswith("PeptideEvidence"):
            # Back to 0, so no cvParam gets written
            stage = 0
        elif entry_tag.endswith("ModificationParams"):
            stage = 0
        elif stage == 1:
            sequence, modifications, peptide_lookup = get_peptide_lookup(
                entry=entry,
                entry_tag=entry_tag,
                sequence=sequence,
                modifications=modifications,
                peptide_lookup=peptide_lookup,
            )
        elif stage == 2:
            modification_mass_map, mod_name, fixed_mods = get_modification_mass_map(
                entry=entry,
                entry_tag=entry_tag,
                modification_mass_map=modification_mass_map,
                mod_name=mod_name,
                fixed_mods=fixed_mods,
            )
        elif stage == 3:
            spec_results, spec_ident_items, spec_records = get_spec_records(
                entry=entry,
                entry_tag=entry_tag,
                spec_results=spec_results,
                spec_ident_items=spec_ident_items,
                spec_records=spec_records,
                mapping_dict=mapping_dict,
            )
        elif entry_tag.endswith("DBSequence"):
            stage = 1
        elif entry_tag.endswith("AdditionalSearchParams"):
            stage = 2
        elif entry_tag.endswith("Inputs"):
            stage = 3
        elif entry_tag.endswith("AnalysisSoftware"):
            version = "comet_" + "_".join(re.findall("[0-9]+", entry.attrib["version"]))
        elif entry_tag.endswith("cvList"):
            if entry_tag != "{http://psidev.info/psi/pi/mzIdentML/1.2}cvList":
                warnings.warn(
                    "Wrong mzIdentML format - Parser might not operate correctly!"
                )
        entry.clear()
        if entry_tag.endswith("SpectrumIdentificationList"):
            return (
                version,
                peptide_lookup,
                spec_records,
                modification_mass_map,
                fixed_mods,
            )


def get_peptide_lookup(entry, entry_tag, sequence, modifications, peptide_lookup):
    """Take one entry at a time to return peptides with their sequences and modifications.

    Args:
        entry (element) : current xml element
        entry_tag (element.tag): was assigned earlier
        sequence (dict): variable from get_xml_data that gets updated
        modifications (str): variable from get_xml_data that gets appended
        peptide_lookup (dict): variable from get_xml_data that gets updated

    Returns:
        sequence (dict): temporary assigned
        modifications (str): temporary assigned
        peptide_lookup (dict): updated dict with one more peptide
    """
    if entry_tag.endswith("PeptideSequence"):
        sequence = {"sequence": entry.text}
    elif entry_tag.endswith("Modification"):
        mod = {"monoisotopicMassDelta": entry.attrib["monoisotopicMassDelta"]}
        mod.update({"location": entry.attrib["location"]})
        modifications.append(mod)
    elif entry_tag.endswith("Peptide"):
        if len(modifications) > 0:
            peptide_lookup[entry.attrib["id"]] = {"modifications": modifications}
        else:
            peptide_lookup[entry.attrib["id"]] = {"modifications": ""}
        peptide_lookup[entry.attrib["id"]].update(sequence)
        modifications = []
    return sequence, modifications, peptide_lookup


def get_modification_mass_map(
    entry, entry_tag, modification_mass_map, mod_name, fixed_mods
):
    """Take one entry at a time to return Modification name with massDelta. Also check if Modification is fixed.

    Args:
        entry (element) : current xml element
        entry_tag (element.tag): was assigned earlier
        modification_mass_map (dict): contains all mods this far
        mod_name (str): contains the name of the Modification
        fixed_mods (dict): contains all fixed mods this far

    Returns:
        modification_mass_map (dict): contains one more modification
        mod_name (str): contains the name of the Modification
        fixed_mods (dict): contains mods where fixedMod = true
    """
    if entry_tag.endswith("cvParam"):
        mod_name = entry.attrib["name"]
    elif entry_tag.endswith("SearchModification"):
        modification_mass_map[entry.attrib["massDelta"]] = mod_name
        if entry.attrib["fixedMod"] == "true":
            fixed_mods.update({entry.attrib["residues"]: mod_name})
    return modification_mass_map, mod_name, fixed_mods


def get_spec_records(
    entry, entry_tag, spec_results, spec_ident_items, spec_records, mapping_dict
):
    """Take one entry at a time to return peptides with their sequences and modifications.

    Args:
        entry (element) : current xml element
        entry_tag (element.tag): was assigned earlier
        spec_results (dict): contains temporary spec information
        spec_ident_items (list): contains all spectrum_identification_items from one spectrum_identification_result
        spec_records (dict): contains all SpectrumIdentificationResult information
        mapping_dict (dict): contains information on which attributes to keep

    Returns:
        spec_results (dict): contains temporary spec information
        spec_ident_items (list): contains all spectrum_identification_items from one spectrum_identification_result
        spec_records (dict): updated with more specs
    """
    if entry_tag.endswith("cvParam"):
        if entry.attrib["name"] in mapping_dict:
            spec_results.update(
                {mapping_dict[entry.attrib["name"]]: entry.attrib["value"]}
            )
    elif entry_tag.endswith("SpectrumIdentificationItem"):
        for attribute in list(entry.attrib):
            if attribute in mapping_dict.keys():
                spec_results.update({mapping_dict[attribute]: entry.attrib[attribute]})
        spec_ident_items.append(spec_results)
        spec_results = {}
    elif entry_tag.endswith("SpectrumIdentificationResult"):
        for spec_item in spec_ident_items:
            spec_item.update(spec_results)
            spec_item.update(
                {"spectrum_id": entry.attrib["spectrumID"].lstrip("scan=")}
            )
            spec_records.append(spec_item)
        spec_results = {}
        spec_ident_items = []
    return spec_results, spec_ident_items, spec_records


class Comet_2020_01_4_Parser(IdentBaseParser):
    """File parser for Comet."""

    def __init__(self, *args, **kwargs):
        """Initialize parser.

        Reads in data file and provides mappings.
        """
        super().__init__(*args, **kwargs)
        self.style = "comet_style_1"
        self.mapping_dict = {
            v: k
            for k, v in self.param_mapper.get_default_params(style=self.style)[
                "header_translations"
            ]["translated_value"].items()
        }

    @classmethod
    def check_parser_compatibility(cls, file):
        """Assert compatibility between file and parser.

        Args:
            file (str): path to input file

        Returns:
            bool: True if parser and file are compatible

        """
        is_mzid = file.as_posix().endswith(".mzid")

        with open(file.as_posix()) as f:
            try:
                head = "".join([next(f) for _ in range(10)])
            except StopIteration:
                head = ""
        contains_engine = "Comet" in head
        return is_mzid and contains_engine

    def map_mods_and_sequences(self, fixed_mods, modification_mass_map, peptide_lookup):
        """Create modification column in df.

        Map correct modifications to corresponding sequence.

        Args:
            fixed_mods (dict): name and symbol of all fixed mods
            modification_mass_map (dict): name and corresponding mass of mods
            peptide_lookup (dict): sequences, with mods from input_file

        Returns:
            self.df (pd.Dataframe): new filled modifications column
        """
        if len(fixed_mods) > 0:
            fixed_mod_strings = []
            for fm_res, fm_name in fixed_mods.items():
                fixed_mod_strings.append(
                    self.df["sequence"]
                    .str.split(fm_res)
                    .apply(
                        lambda l: ";".join(
                            [
                                fm_name + ":" + ind
                                for ind in (
                                    np.cumsum(list(map(len, l[:-1]))) + range(1, len(l))
                                ).astype(str)
                            ]
                        )
                    )
                )

            fixed_mod_strings = (
                pd.concat(fixed_mod_strings, axis=1)
                .agg(";".join, axis=1)
                .str.rstrip(";")
            )
        lookup = {}
        for pep_sequence, attribs in peptide_lookup.items():
            lookup[pep_sequence] = {
                "modifications": [],
                "sequence": attribs["sequence"],
            }
            if len(attribs["modifications"]) != 0:
                for mod in attribs["modifications"]:
                    monoisotopicMassDelta = mod["monoisotopicMassDelta"]
                    location = mod["location"]
                    lookup[pep_sequence]["modifications"].append(
                        f"{modification_mass_map[monoisotopicMassDelta]}:{location}"
                    )
            lookup[pep_sequence]["modifications"] = ";".join(
                lookup[pep_sequence]["modifications"]
            )

        seq_mods = pd.DataFrame(self.df["sequence"].map(lookup).to_list())
        self.df.loc[:, "modifications"] = (
            seq_mods["modifications"].str.cat(fixed_mod_strings, sep=";").str.strip(";")
        )
        self.df.loc[:, "sequence"] = seq_mods["sequence"]

    def unify(self):
        """
        Primary method to read and unify engine output.

        Returns:
            self.df (pd.DataFrame): unified dataframe
        """
        (
            version,
            peptide_lookup,
            spec_records,
            modification_mass_map,
            fixed_mods,
        ) = get_xml_data(self.input_file, self.mapping_dict)
        self.df = pd.DataFrame(spec_records)
        self.df["search_engine"] = version
        self.map_mods_and_sequences(fixed_mods, modification_mass_map, peptide_lookup)
        self.process_unify_style()

        return self.df
