"""Engine parser."""
import xml.etree.ElementTree as etree

import pandas as pd
import regex as re
from loguru import logger

from pyiohat.parsers.ident_base_parser import IdentBaseParser


class MSGFPlus_2021_03_22_Parser(IdentBaseParser):
    """File parser for MSGF+."""

    def __init__(self, *args, **kwargs):
        """Initialize parser.

        Reads in data file and provides mappings.
        """
        super().__init__(*args, **kwargs)
        self.version = None
        self.peptide_lookup = None
        self.spec_records = None
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
            file (path object): path to input file

        Returns:
            bool: True if parser and file are compatible
        """
        is_mzid = file.name.endswith(".mzid")

        with open(file) as f:
            try:
                head = "".join([next(f) for _ in range(20)])
            except StopIteration:
                head = ""
        contains_engine = "MS-GF+" in head
        return is_mzid and contains_engine

    def get_version(self):
        """Retrieve version from xml.

        Returns:
            version (str): file version
        """
        version = ""
        for event, entry in etree.iterparse(self.input_file):
            entry_tag = entry.tag

            if entry_tag.endswith("cvList"):
                if entry_tag != "{http://psidev.info/psi/pi/mzIdentML/1.1}cvList":
                    logger.warning(
                        f"{entry_tag}: Wrong mzIdentML version - Parser made for version 1.1.0!"
                    )
            elif entry_tag.endswith("AnalysisSoftware"):
                version = "msgfplus_" + "_".join(
                    re.findall("[0-9]+", entry.attrib["version"])
                )
                break
            entry.clear()
        return version

    def get_peptide_lookup(self):
        """Retrieve peptide ids with their corresponding sequences and modifications from xml.

        Returns:
            peptide_lookup (dict): peptide id with corresponding modifications and sequence
        """
        peptide_lookup = {}

        cv_param_modifications = ""
        peptide_information = False

        for event, entry in etree.iterparse(self.input_file):
            entry_tag = entry.tag

            if entry_tag.endswith("PeptideSequence"):
                peptide_information = True
            if peptide_information is True:
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
                    peptide_lookup[entry.attrib["id"]] = sequence
                    peptide_lookup[entry.attrib["id"]][
                        "modifications"
                    ] = cv_param_modifications.rstrip(";")
                    cv_param_modifications = ""
                elif entry_tag.endswith("PeptideEvidence"):
                    break
            entry.clear()
        return peptide_lookup

    def get_spec_records(self):
        """Retrieve specs from file.

        Returns:
            spec_records (list): information on PSMs
        """
        spec_records = []

        spec_results = {}
        spec_ident_items = []
        spec_information = False

        for event, entry in etree.iterparse(self.input_file):
            entry_tag = entry.tag

            if entry_tag.endswith("FragmentationTable"):
                spec_information = True
            elif spec_information is True:
                if entry_tag.endswith("cvParam") or entry_tag.endswith("userParam"):
                    if entry.attrib["name"] in self.mapping_dict:
                        _key = self.mapping_dict[entry.attrib["name"]]
                        spec_results[_key] = entry.attrib["value"]
                elif entry_tag.endswith("SpectrumIdentificationItem"):
                    for attribute in list(entry.attrib):
                        if attribute in self.mapping_dict.keys():
                            if attribute == "peptide_ref":
                                psm_id = entry.attrib[attribute]
                                spec_results["sequence"] = self.peptide_lookup[psm_id][
                                    "sequence"
                                ]
                                spec_results["modifications"] = self.peptide_lookup[
                                    psm_id
                                ]["modifications"]
                            else:
                                _key = self.mapping_dict[attribute]
                                spec_results[_key] = entry.attrib[attribute]
                    # multiple SpectrumIdentificationItems possible, therefore create a list and reset spec_results
                    spec_ident_items.append(spec_results)
                    spec_results = {}
                elif entry_tag.endswith("SpectrumIdentificationResult"):
                    for spec_item in spec_ident_items:
                        spec_item.update(spec_results)
                        spec_records.append(spec_item)
                    spec_results = {}
                    spec_ident_items = []
                elif entry_tag.endswith("SpectrumIdentificationList"):
                    break
            entry.clear()
        return spec_records

    def unify(self):
        """
        Primary method to read and unify engine output.

        Returns:
            self.df (pd.DataFrame): unified dataframe
        """
        self.version = self.get_version()
        self.peptide_lookup = self.get_peptide_lookup()
        self.spec_records = self.get_spec_records()
        self.df = pd.DataFrame(self.spec_records)
        self.df["search_engine"] = self.version
        self.process_unify_style()

        return self.df
