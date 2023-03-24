"""Quant parser."""
import re
from pathlib import Path

import pandas as pd
from chemical_composition import ChemicalComposition

from pyiohat.parsers.quant_base_parser import QuantBaseParser
from pyiohat.parsers.misc import get_compositions_and_monoisotopic_masses


class TMTQuantParser(QuantBaseParser):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.style = "cz_tmt_quant_style_1"
        self.df = pd.read_csv(self.input_file)

        self.mapping_dict = {
            v: k
            for k, v in self.param_mapper.get_default_params(style=self.style)[
                "header_translations"
            ]["translated_value"].items()
        }
        self.df.rename(columns=self.mapping_dict, inplace=True)

    @classmethod
    def check_parser_compatibility(cls, file):
        is_csv = file.as_posix().endswith(".csv")
        with open(file.as_posix(), "r") as fin:
            try:
                header = "".join([next(fin) for _ in range(1)])
            except StopIteration:
                header = ""
        header = set(header.rstrip("\n").split(","))

        ref_columns = {
            "filename",
            "mz",
            "raw_quant_intensity",
            "retention_time_seconds",
            "label",
            "mz_delta",
            "iso_mz",
            "ppm",
            "spectrum_id",
            "raw_quant_area",
            "isolabel_id",
            "original_quant_value",
            "s2i",
            "quant_value",
        }
        columns_match = len(ref_columns.difference(header)) == 0
        return columns_match and is_csv

    def unify(self):
        """
        Primary method to read and unify engine output.

        Returns:
            self.df (pd.DataFrame): unified dataframe
        """
        # base model
        self.df["charge"] = -1
        self.df["accuracy_ppm"] = -1
        self.df["accuracy_ppm_C12"] = -1
        self.df["chemical_composition"] = ""
        # quant model
        self.df["trivial_name"] = ""
        self.df["quant_score"] = -1
        self.df["quant_group"] = ""
        self.df["ccs"] = -1
        self.df["ident_reference"] = -1
        self.df["fwhm"] = -1
        self.df["quant_engine"] = "TMTQuant_1_0_0"
        self.process_unify_style()
        # required columns for recalculation not in dataframe
        self.df["accuracy_ppm"] = -1
        return self.df
