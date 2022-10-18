"""Quant base parser class."""

from chemical_composition import ChemicalComposition

from pyprotista.parsers.base_parser import BaseParser


class QuantBaseParser(BaseParser):
    """Base class of all quant parsers."""

    def __init__(self, *args, **kwargs):
        """Initialize parser.

        Reads in data file and provides mappings.
        """
        super().__init__(*args, **kwargs)
        self.cc = ChemicalComposition()
        self.required_headers = self._load_model("quant_parser_model.json")

    def process_unify_style(self):
        """Apply sanitizing methods."""
        # TODO: this should be called here too
        # self.sanitize()
