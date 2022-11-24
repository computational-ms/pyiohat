"""Engine parser."""
import pandas as pd
import regex as re
import xml.etree.ElementTree as etree
from pyprotista.parsers.ident_base_parser import IdentBaseParser
from loguru import logger


def get_xml_data(file, mapping_dict):
    """Iterate over file to get information on specs.

    Args:
        file (str): path to input file
        mapping_dict (dict): mapping of engine level column names to pyprotista unified column names

    Returns:
        version (str): file version
        peptide_lookup (dict): peptide id with corresponding modifications and sequence
        spec_records (list): spectrum information
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

    for event, entry in etree.iterparse(file):
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
                logger.warning(
                    "Wrong mzIdentML format - Parser might not operate correctly!"
                )
        entry.clear()
        if entry_tag.endswith("SpectrumIdentificationList"):
            return version, peptide_lookup, spec_records


def get_peptide_lookup(
    entry, entry_tag, sequence, cv_param_modifications, peptide_lookup
):
    """Take one entry at a time to get peptide ids with their corresponding sequences and modifications.

    Args:
        entry (element): xml element
        entry_tag (element.tag): xml tag
        sequence (dict): current sequence
        cv_param_modifications (str): modifications of corresponding sequence
        peptide_lookup (dict): peptide id with corresponding modifications and sequence

    Returns:
        sequence (dict): current sequence
        cv_paparam_modifications (str): temporary assigned
        peptide_lookup (dict): peptide id with corresponding modifications and sequence
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
    """Take one entry at a time to get single specs.

    Args:
        entry (element) : xml element
        entry_tag (str): xml tag
        spec_ident_items (list): potentially multiple PSMs
        spec_results (dict): single PSM
        spec_records (list): information on PSMs
        mapping_dict (dict): mapping of engine level column names to pyprotista unified column names

    Returns:
        spec_results (dict): single PSM
        spec_ident_items (list): potentially multiple PSMs
        spec_records (list): information on PSMs
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
