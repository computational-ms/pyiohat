"""Quant base parser class."""

from pyiohat.parsers.base_parser import BaseParser
import pandas as pd


class QuantBaseParser(BaseParser):
    """Base class of all quant parsers."""

    def __init__(self, *args, **kwargs):
        """Initialize parser.

        Reads in data file and provides mappings.
        """
        super().__init__(*args, **kwargs)
        self.required_headers = self._load_model("quant_parser_model.json")
        self.col_order = pd.Series(self.required_headers.keys())
        self.rt_lookup = self._read_meta_info_lookup_file()

    def process_unify_style(self):
        """Apply sanitizing methods."""
        # TODO: this should be called here too
        self.calculate_accuracy()
        self.sanitize()

    def calculate_accuracy(self):
        accuracy_mz = self.df.loc[
            self.df["reported_mz"] != "-", "theoretical_mz"
        ].astype(float) - self.df.loc[
            self.df["reported_mz"] != "-", "reported_mz"
        ].astype(
            float
        )
        accuracy_ppm = (
            (accuracy_mz)
            / self.df.loc[self.df["reported_mz"] != "-", "reported_mz"].astype(float)
        ) * 1e6
        self.df.loc[self.df["reported_mz"] != "-", "accuracy_ppm"] = accuracy_ppm
        self.df.loc[self.df["reported_mz"] == "-", "accuracy_ppm"] = -1
        self.df.loc[self.df["reported_mz"] == "-", "reported_ppm"] = -1

        self.df.loc[self.df["reported_mz"] != "-", "accuracy_mz"] = accuracy_mz
        self.df.loc[self.df["reported_mz"] == "-", "accuracy_mz"] = -1
        self.df.loc[self.df["reported_mz"] == "-", "reported_mz"] = -1
