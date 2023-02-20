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
        # self.sanitize()
        pass
