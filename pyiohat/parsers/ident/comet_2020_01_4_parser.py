"""Engine parser."""
import xml.etree.ElementTree as etree

import numpy as np
import pandas as pd
import regex as re
from loguru import logger

from pyiohat.parsers.ident_base_parser import IdentBaseParser


class Comet_2020_01_4_Parser(IdentBaseParser):
    """File parser for Comet."""

    def __init__(self, *args, **kwargs):
        """Initialize parser.

        Reads in data file and provides mappings.
        """
        super().__init__(*args, **kwargs)
        self.version = None
        self.fixed_mods = None
        self.mod_mass_map = None
        self.peptide_lookup = None
        self.spec_records = None
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
            file (path object): path to input file

        Returns:
            bool: True if parser and file are compatible
        """
        is_mzid = file.name.endswith(".mzid")

        with open(file) as f:
            try:
                head = "".join([next(f) for _ in range(10)])
            except StopIteration:
                head = ""
        contains_engine = "Comet" in head
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
                if entry_tag != "{http://psidev.info/psi/pi/mzIdentML/1.2}cvList":
                    logger.warning(
                        f"{entry_tag}: Wrong mzIdentML version - Parser made for version 1.2!"
                    )
            elif entry_tag.endswith("AnalysisSoftware"):
                version = "comet_" + "_".join(
                    re.findall("[0-9]+", entry.attrib["version"])
                )
                break
            entry.clear()
        return version

    def map_mod_mass(self):
        """Retrieve information on modifications from xml.

        Returns:
            fixed_mods (dict): residues and corresponding names of fixed modifications
            mod_mass_map (dict): mapping of masses to modification names
        """
        fixed_mods = {}
        mod_mass_map = {}

        mod_name = ""
        modification_information = False

        for event, entry in etree.iterparse(self.input_file):
            entry_tag = entry.tag

            if entry_tag.endswith("AdditionalSearchParams"):
                modification_information = True
            elif modification_information is True:
                if entry_tag.endswith("cvParam"):
                    mod_name = entry.attrib["name"]
                elif entry_tag.endswith("SearchModification"):
                    if mod_name == "unknown modification":
                        potential_mod = self.mod_mapper.mass_to_names(
                            float(entry.attrib["massDelta"]), decimals=4
                        )
                        if len(potential_mod) == 0:
                            logger.error(
                                f"Cannot map modification with mass {entry.attrib['massDelta']}."
                            )
                            raise ValueError
                        else:
                            mod_name = potential_mod[0]
                    mod_mass_map[entry.attrib["massDelta"]] = mod_name
                    if entry.attrib["fixedMod"] == "true":
                        residue = entry.attrib["residues"]
                        fixed_mods[residue] = mod_name
                elif entry_tag.endswith("ModificationParams"):
                    break
            entry.clear()
        return fixed_mods, mod_mass_map

    def get_peptide_lookup(self):
        """Retrieve peptide ids with their corresponding sequences and modifications from xml.

        Returns:
            peptide_lookup (dict): peptide id with corresponding modifications and sequence
        """
        peptide_lookup = {}

        modifications = []
        sequence = {}
        peptide_information = False

        for event, entry in etree.iterparse(self.input_file):
            entry_tag = entry.tag

            if entry_tag.endswith("PeptideSequence"):
                peptide_information = True
            if peptide_information is True:
                if entry_tag.endswith("PeptideSequence"):
                    sequence = entry.text
                    if len(self.fixed_mods) > 0:
                        sequence_mod_map = self.map_mods_sequences(sequence)
                        if sequence_mod_map != "":
                            modifications.append(sequence_mod_map)
                elif entry_tag.endswith("Modification"):
                    mass = entry.attrib["monoisotopicMassDelta"]
                    location = entry.attrib["location"]
                    mass_name = self.mod_mass_map[mass]
                    modifications.append(mass_name + ":" + location)
                elif entry_tag.endswith("Peptide"):
                    peptide_lookup[entry.attrib["id"]] = (
                        ";".join(modifications),
                        sequence,
                    )
                    modifications = []
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

            if entry_tag.endswith("PeptideEvidenceRef"):
                spec_information = True
            if spec_information is True:
                if entry_tag.endswith("cvParam"):
                    if entry.attrib["name"] in self.mapping_dict:
                        _key = self.mapping_dict[entry.attrib["name"]]
                        spec_results[_key] = entry.attrib["value"]
                elif entry_tag.endswith("SpectrumIdentificationItem"):
                    for attribute in list(entry.attrib):
                        if attribute in self.mapping_dict.keys():
                            _key = self.mapping_dict[attribute]
                            spec_results[_key] = entry.attrib[attribute]
                    mods, sequence = self.peptide_lookup[spec_results["sequence"]]
                    spec_results["modifications"] = mods
                    spec_results["sequence"] = sequence
                    spec_ident_items.append(spec_results)
                    spec_results = {}
                elif entry_tag.endswith("SpectrumIdentificationResult"):
                    for spec_item in spec_ident_items:
                        spec_item.update(spec_results)
                        spec_item["spectrum_id"] = entry.attrib["spectrumID"].lstrip(
                            "scan="
                        )
                        spec_records.append(spec_item)
                    spec_results = {}
                    spec_ident_items = []
                elif entry_tag.endswith("SpectrumIdentificationList"):
                    break
            entry.clear()
        return spec_records

    def map_mods_sequences(self, sequence):
        """Map fixed_mods and sequence to get corresponding modifications.

        fixed_mods needs at least one entry!

        Args:
            sequence (str): peptide sequence

        Returns:
            fixed_mod_strings (str): modifications of corresponding sequence
        """
        fixed_mod_strings = []

        for fm_res, fm_name in self.fixed_mods.items():
            fixed_mod_strings.append(
                pd.Series(sequence)
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
            pd.concat(fixed_mod_strings, axis=1).agg(";".join, axis=1).str.rstrip(";")
        )
        fixed_mod_strings = fixed_mod_strings[0]
        return fixed_mod_strings

    def unify(self):
        """
        Primary method to read and unify engine output.

        Returns:
            self.df (pd.DataFrame): unified dataframe
        """
        self.version = self.get_version()
        self.fixed_mods, self.mod_mass_map = self.map_mod_mass()
        self.peptide_lookup = self.get_peptide_lookup()
        self.spec_records = self.get_spec_records()
        self.df = pd.DataFrame(self.spec_records)
        self.df["search_engine"] = self.version
        self.process_unify_style()

        return self.df
