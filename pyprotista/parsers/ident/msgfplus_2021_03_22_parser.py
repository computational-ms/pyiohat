"""Engine parser."""
import pandas as pd
import regex as re
import xml.etree.cElementTree as etree
from pyprotista.parsers.ident_base_parser import IdentBaseParser

element_tag_prefix = "{http://psidev.info/psi/pi/mzIdentML/1.1}"


def _iterator_xml(xml_file, mapping_dict):
    version = ""
    stage = 0
    cv_param_modifications = ""
    sequence = {}
    peptide_lookup = {}
    spec_results = {}
    spec_ident_items = []
    spec_records = []

    for event, entry in etree.iterparse(xml_file):
        entry_tag = entry.tag

        if entry_tag == (f"{element_tag_prefix}PeptideEvidence"):
            stage = 0
            continue
            # Back to 0, so no cvParam gets written

        elif stage == 1:
            sequence, cv_param_modifications, peptide_lookup = _peptide_lookup(
                entry, entry_tag, sequence, cv_param_modifications, peptide_lookup
            )
        elif stage == 2:
            spec_results, spec_ident_items, spec_records = _spec_records(
                entry,
                entry_tag,
                spec_results,
                spec_ident_items,
                spec_records,
                mapping_dict,
            )

        elif entry_tag == (f"{element_tag_prefix}DBSequence"):
            stage = 1
            continue
            # Short before peptide_lookup begins
        elif entry_tag == (f"{element_tag_prefix}FragmentationTable"):
            stage = 2
            continue
            # Short before spec_record begins

        elif entry_tag == (f"{element_tag_prefix}AnalysisSoftware"):
            version = "msgfplus_" + "_".join(
                re.findall(r"([/d]*\d+)", entry.attrib["version"])
            )

        entry.clear()
        if entry_tag == (f"{element_tag_prefix}SpectrumIdentificationList"):
            return version, peptide_lookup, spec_records


def _peptide_lookup(entry, entry_tag, sequence, cv_param_modifications, peptide_lookup):
    if entry_tag == (f"{element_tag_prefix}PeptideSequence"):
        sequence = {"sequence": entry.text}
    elif entry_tag == (f"{element_tag_prefix}cvParam"):
        if entry.attrib["name"] == "unknown modification":
            cv_param_modifications += entry.attrib["value"] + ":"
        else:
            cv_param_modifications += entry.attrib["name"] + ":"
    elif entry_tag == (f"{element_tag_prefix}Modification"):
        cv_param_modifications += entry.attrib["location"] + ";"
    elif entry_tag == (f"{element_tag_prefix}Peptide"):
        peptide_lookup[entry.attrib["id"]] = {
            "modifications": cv_param_modifications.rstrip(";")
        }
        peptide_lookup[entry.attrib["id"]].update(sequence)
        cv_param_modifications = ""
        sequence = ""
    return sequence, cv_param_modifications, peptide_lookup


def _spec_records(
    entry, entry_tag, spec_results, spec_ident_items, spec_records, mapping_dict
):
    if entry_tag == (f"{element_tag_prefix}cvParam") or entry_tag == (
        f"{element_tag_prefix}userParam"
    ):
        if entry.attrib["name"] in mapping_dict:
            spec_results.update(
                {mapping_dict[entry.attrib["name"]]: entry.attrib["value"]}
            )
    elif entry_tag == (f"{element_tag_prefix}SpectrumIdentificationItem"):
        for i in list(entry.attrib):
            if i in mapping_dict:
                spec_results.update({mapping_dict[i]: entry.attrib[i]})
        spec_ident_items.append(spec_results)
        spec_results = {}
    elif entry_tag == (f"{element_tag_prefix}SpectrumIdentificationResult"):
        for i in spec_ident_items:
            i.update(spec_results)
            spec_records.append(i)
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
        self.reference_dict.update({k: None for k in self.mapping_dict.values()})

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
        version, peptide_lookup, spec_records = _iterator_xml(
            self.input_file, self.mapping_dict
        )
        self.df = pd.DataFrame(spec_records)
        seq_mods = pd.DataFrame(self.df["sequence"].map(peptide_lookup).to_list())
        self.df.loc[:, seq_mods.columns] = seq_mods
        self.df["search_engine"] = version
        self.process_unify_style()

        return self.df
