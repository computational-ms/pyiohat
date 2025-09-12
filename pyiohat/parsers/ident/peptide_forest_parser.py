"""Engine parser."""

import itertools

import pandas as pd
import regex as re
from loguru import logger
from importlib import import_module
from pyiohat.parsers.ident_base_parser import IdentBaseParser
from pyiohat.parsers.base_parser import BaseParser
from pprint import pprint
from itertools import combinations
from chemical_composition import chemical_composition_kb
from pathlib import Path


class PeptideForest_Parser(IdentBaseParser):
    """File parser for PTMShepherd"""

    def __init__(self, *args, **kwargs):
        """Initialize parser.

        Reads in data file and provides mappings.
        """
        super().__init__(*args, **kwargs)
        self.style = "peptide_forest_style_1"

        self.df = pd.read_csv(self.input_file, delimiter=",")
        self.df.dropna(axis=1, how="all", inplace=True)
        original_columns = set(self.df.columns)
        untouched_columns = [
            "raw_data_location",
            "spectrum_id",
            "sequence",
            "modifications",
            "is_decoy",
            "protein_id",
            "charge",
            "accuracy_ppm",
            "accuracy_ppm_C12",
            "enzc",
            "enzn",
            "exp_mass",
            "is_immutable",
            "mass_delta",
            "missed_cleavages",
            "retention_time_seconds",
            "ucalc_mass",
            "pep_len",
            "count_prot",
        ]
        unmapped_columns = original_columns.difference(untouched_columns)
        prefix_mapping_dict = {col: f"peptide_forest:{col}" for col in unmapped_columns}
        self.df.rename(columns=prefix_mapping_dict, inplace=True)
        self.df.rename(columns=self.mapping_dict, inplace=True)
        self.df.columns = self.df.columns.str.lstrip(" ")
        self.reference_dict.update({k: None for k in self.mapping_dict.values()})
        self.metadata = {
            "File Origin": "PeptideForest",
            "Version": ["3.1"],
            "bigger_scores_better": True,
            "validation_score_field": "peptide_forest:q-value_peptide_forest",
            "Parser": "pyiohat/parsers/ident/peptide_forest_parser.py",
        }

    @classmethod
    def check_parser_compatibility(cls, file):
        """Assert compatibility between file and parser.

        Args:
            file (str): path to input file

        Returns:
            bool: True if parser and file are compatible

        """
        is_csv = file.as_posix().endswith(".csv")
        with open(file.as_posix()) as f:
            try:
                head = "".join([next(f) for _ in range(1)])
            except StopIteration:
                head = ""
        head = set(head.rstrip("\n").split("\t"))
        ref_columns = {
            "raw_data_location",
            "spectrum_id",
            "sequence",
            "modifications",
            "is_decoy",
            "protein_id",
            "charge",
            "accuracy_ppm",
            "accuracy_ppm_C12",
            "enzc",
            "enzn",
            "exp_mass",
            "is_immutable",
            "mass_delta",
            "missed_cleavages",
            "retention_time_seconds",
            "ucalc_mass",
            "pep_len",
            "count_prot",
            "score_processed_peptide_forest",
            "q-value_peptide_forest",
            "top_target_peptide_forest",
            "rank_peptide_forest",
        }
        columns_match = len(ref_columns.difference(head)) == 0
        return is_csv and columns_match

    def unify(self):
        """
        Primary method to read and unify engine output.

        Returns:
            self.df (pd.DataFrame): unified dataframe
        """
        self.df["validation_engine"] = "peptide_forest_3"
        self.process_unify_style()

        return self.df
