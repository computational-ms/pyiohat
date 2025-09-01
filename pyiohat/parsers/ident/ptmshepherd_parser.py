"""Engine parser."""

import itertools

import pandas as pd
import regex as re
from loguru import logger
from importlib import import_module
from pyiohat.parsers.ident_base_parser import IdentBaseParser
from pprint import pprint
from itertools import combinations
from chemical_composition import chemical_composition_kb
from pathlib import Path


class PTMShepherd_Parser(IdentBaseParser):
    """File parser for PTMShepherd"""

    def __init__(self, *args, **kwargs):
        """Initialize parser.

        Reads in data file and provides mappings.
        """
        super().__init__(*args, **kwargs)
        self.style = "ptmshepherd_style_1"
        # 15N handling missing for now
        if self.params.get("label", "") == "15N":
            raise NotImplementedError

        self.df = pd.read_csv(self.input_file, delimiter="\t")
        self.df.dropna(axis=1, how="all", inplace=True)
        self.mapping_dict = {
            "Spectrum File": "raw_data_location",
            "Charge": "charge",
            "Retention": "retention_time_seconds",
            "Spectrum": "spectrum_title",
            "Peptide": "sequence",
            "Assigned Modifications": "modifications",
            "Protein ID": "protein_id",
            "Observed M/Z": "exp_mz",
            "Observed Mass": "exp_mass",
            "Calculated M/Z": "ucalc_mz",
            "Calculated Peptide Mass": "ucalc_mass",
            "Delta Mass": "mass_delta",
            "MSFragger Localization": "Localization",
        }
        self.df.rename(columns=self.mapping_dict, inplace=True)
        self.df.columns = self.df.columns.str.lstrip(" ")
        if not "modifications" in self.df.columns:
            self.df["modifications"] = ""
        self.reference_dict.update({k: None for k in self.mapping_dict.values()})
        self.metadata = self._get_metadata()


    def _get_metadata(self):
        
        metadata = {
            "File Origin": "PTMShepherd",
            "Version": [2.0.5],
            "Parser": "pyiohat/parsers/ident/ptmshepherd_parser.py",
        }

        if not "search_engine" in self.df.columns:
            # This means the parser is running directly after search + ptmshepherd without pyiohat inbetween
            # search engine must have been msfragger, write bigger_score_better and validation_score_field accordinglly
            metadata.append("validation_score_field": "ptmshepherd:hyperscore", "bigger_scores_better": True)

        else:
            # This means the parser is running after search + pyiohat + ptmshepherd
            # In this case bigger_score_better and validation_score_field needs to be retrived from the parser of the relavent search engine
            parsers_dict = {
            "xtandem_":"xtandem_alanine",
            "omssa_2_1_9":"omssa_2_1_9_parser",
            "msgfplus_":"msgfplus_2021_03_22_parser",
            "msfragger_4_2": "msfragger_4_parser",
            "msfragger_3_0":"msfragger_3_parser",
            "msamanda_2_0_0_17442":"msamanda_2_parser",
            "mascot_":"mascot_2_6_2_parser",
            "comet_":"comet_2020_01_4_parser",
            }
            search_engine = self.df["search_engine"][1]
            for k, v in parsers_dict.items():
                if k in search_engine:
                    parser_name = v
            original_parser_module = f"{__package__}.{module_name}"
            import_module(original_parser_module)
            parser_classes = []
            for cat in BaseParser.__subclasses__():
                parser_classes.extend(
                    [c for c in cat.__subclasses__() if c.__module__ == original_parser_module]
                )
            ParserClass = parser_classes[0]
            parser_instance = ParserClass(
                input_file=self.input_file,
                params=self.params,
                immutable_peptides=self.immutable_peptides,
            )
            
            original_metadata = parser_instance.metadata

            metadata.append("validation_score_field": original_metadata["validation_score_field"], "bigger_scores_better": original_metadata["bigger_scores_better"])

        return metadata

 
    @classmethod
    def check_parser_compatibility(cls, file):
        """Assert compatibility between file and parser.

        Args:
            file (str): path to input file

        Returns:
            bool: True if parser and file are compatible

        """
        is_tsv = file.as_posix().endswith(".tsv")
        with open(file.as_posix()) as f:
            try:
                head = "".join([next(f) for _ in range(1)])
            except StopIteration:
                head = ""
        head = set(head.rstrip("\n").split("\t"))
        ref_columns = {
            "Spectrum File",
            "Charge",
            "Retention",
            "Spectrum",
            "Peptide",
            "Assigned Modifications",
            "Protein ID",
            "Observed M/Z",
            "Observed Mass",
            "Calculated M/Z",
            "Calculated Peptide Mass",
            "Delta Mass",
            "MSFragger Localization",
            "Modified Peptide",
            "Observed Modifications",            
        }
        columns_match = len(ref_columns.difference(head)) == 0
        return is_tsv and columns_match

    def unify(self):
        """
        Primary method to read and unify engine output.

        Returns:
            self.df (pd.DataFrame): unified dataframe
        """
        self.df["glycan_annotation_engine"] = "ptmshepherd"
        self.process_unify_style()

        return self.df
