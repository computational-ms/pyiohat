"""Engine parser."""
import pandas as pd
import regex as re
import xml.etree.ElementTree as etree
from pyprotista.parsers.ident_base_parser import IdentBaseParser
import warnings


def get_xml_data(xml_file, mapping_dict):
    """For loop over one xml file using xml.etree.ElementTree.iterparse.

    Provide temporary variables for iteration.
    Call functions get_peptide_lookup in stage 1 and get_spec_records in stage 2.

    Args:
        xml_file (mzid): input file
        mapping_dict (dict)

    Returns:
        version (str): contains the version of the mzid
        peptide_lookup (dict): contains all the peptides with their sequences and modifications
        spec_records (list): list of dicts containing single_specs
    """
    version = ""
    peptide_lookup = {}
    spec_records = []

    # temporary variables, get overwritten multiple times during iteration
    stage = 0
    cv_param_modifications = ""
    sequence = {}
    spec_results = {}
    spec_ident_items = []

    for event, entry in etree.iterparse(xml_file):
        entry_tag = entry.tag
        if entry_tag.endswith("PeptideEvidence"):
            # Back to 0, so no cvParam gets written
            stage = 0
        elif stage == 1:
            sequence, cv_param_modifications, peptide_lookup = get_peptide_lookup(
                entry=entry,
                entry_tag=entry_tag,
                sequence=sequence,
                cv_param_modifications=cv_param_modifications,
                peptide_lookup=peptide_lookup,
            )
        elif stage == 2:
            spec_results, spec_ident_items, spec_records = get_spec_records(
                entry=entry,
                entry_tag=entry_tag,
                spec_ident_items=spec_ident_items,
                spec_results=spec_results,
                spec_records=spec_records,
                mapping_dict=mapping_dict,
            )
        elif entry_tag.endswith("DBSequence"):
            stage = 1
        elif entry_tag.endswith("FragmentationTable"):
            stage = 2
        elif entry_tag.endswith("AnalysisSoftware"):
            version = "msgfplus_" + "_".join(
                re.findall("[0-9]+", entry.attrib["version"])
            )
        elif entry_tag.endswith("cvList"):
            if entry_tag != "{http://psidev.info/psi/pi/mzIdentML/1.1}cvList":
                warnings.warn(
                    "Wrong mzIdentML format - Parser might not operate correctly!"
                )
        # RAM saver
        entry.clear()
        if entry_tag.endswith("SpectrumIdentificationList"):
            return version, peptide_lookup, spec_records


def get_peptide_lookup(
    entry, entry_tag, sequence, cv_param_modifications, peptide_lookup
):
    """Take one entry at a time to return peptides with their sequences and modifications.

    Args:
        entry (element) : current xml element
        entry_tag (element.tag): was assigned earlier
        sequence (dict): variable from get_xml_data that gets updated
        cv_param_modifications (str): variable from get_xml_data that gets appended
        peptide_lookup (dict): variable from get_xml_data that gets updated

    Returns:
        sequence (dict): temporary assigned
        cv_paparam_modifications (str): temporary assigned
        peptide_lookup (dict): updated dict with one more peptide
    """
    if entry_tag.endswith("PeptideSequence"):
        sequence = {"sequence": entry.text}
    elif entry_tag.endswith("cvParam"):
        if entry.attrib["name"] == "unknown modification":
            cv_param_modifications += entry.attrib["value"] + ":"
        else:
            cv_param_modifications += entry.attrib["name"] + ":"
    elif entry_tag.endswith("Modification"):
        cv_param_modifications += entry.attrib["location"] + ";"
    elif entry_tag.endswith("Peptide"):
        peptide_lookup[entry.attrib["id"]] = {
            "modifications": cv_param_modifications.rstrip(";")
        }
        peptide_lookup[entry.attrib["id"]].update(sequence)
        cv_param_modifications = ""
    return sequence, cv_param_modifications, peptide_lookup


def get_spec_records(
    entry, entry_tag, spec_ident_items, spec_results, spec_records, mapping_dict
):
    """Take one entry at a time to return peptides with their sequences and modifications.

    Args:
        entry (element) : current xml element
        entry_tag (element.tag): was assigned earlier
        spec_ident_items (list): contains all spectrum_identification_items from one spectrum_identification_result
        spec_results (dict): contains temporary spec information
        spec_records (dict): contains all SpectrumIdentificationResult information
        mapping_dict (dict): contains information on which attributes to keep

    Returns:
        spec_results (dict): contains temporary spec information
        spec_ident_items (list): contains all spectrum_identification_items from one spectrum_identification_result
        spec_records (dict): updated with more specs
    """
    if entry_tag.endswith("cvParam") or entry_tag.endswith("userParam"):
        if entry.attrib["name"] in mapping_dict:
            spec_results.update(
                {mapping_dict[entry.attrib["name"]]: entry.attrib["value"]}
            )
    elif entry_tag.endswith("SpectrumIdentificationItem"):
        for attribute in list(entry.attrib):
            if attribute in mapping_dict.keys():
                spec_results.update({mapping_dict[attribute]: entry.attrib[attribute]})
        # multiple SpectrumIdentificationItems possible, therefore create a list and reset spec_results
        spec_ident_items.append(spec_results)
        spec_results = {}

    elif entry_tag.endswith("SpectrumIdentificationResult"):
        for spec_item in spec_ident_items:
            spec_item.update(spec_results)
            spec_records.append(spec_item)
        spec_results = {}
        spec_ident_items = []
    return spec_results, spec_ident_items, spec_records


class MSGFPlus_2021_03_22_Parser(IdentBaseParser):
    """File parser for MSGF+."""

    def __init__(self, *args, **kwargs):
        """Initialize parser.

        Reads in data file and provides mappings.
        """
        super().__init__(*args, **kwargs)
        self.style = "msgfplus_style_1"
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
        is_mzid = file.name.endswith(".mzid")

        with open(file.as_posix()) as f:
            try:
                head = "".join([next(f) for _ in range(20)])
            except StopIteration:
                head = ""
        contains_engine = "MS-GF+" in head
        return is_mzid and contains_engine

    def unify(self):
        """
        Primary method to read and unify engine output.

        Returns:
            self.df (pd.DataFrame): unified dataframe
        """
        version, peptide_lookup, spec_records = get_xml_data(
            self.input_file, self.mapping_dict
        )
        self.df = pd.DataFrame(spec_records)
        seq_mods = pd.DataFrame(self.df["sequence"].map(peptide_lookup).to_list())
        self.df.loc[:, seq_mods.columns] = seq_mods
        self.df["search_engine"] = version
        self.process_unify_style()

        return self.df
