"""Engine parser."""

import pandas as pd
import regex as re
from loguru import logger

from pyiohat.parsers.ident_base_parser import IdentBaseParser
from pprint import pprint


class GlycoDecipher_1_Parser(IdentBaseParser):
    """File parser for pGlyco 3"""

    def __init__(self, *args, **kwargs):
        """Initialize parser.

        Reads in data file and provides mappings.
        """
        super().__init__(*args, **kwargs)
        self.style = "glyco_decipher_style_1"

        self.df = pd.read_csv(self.input_file, delimiter="\t")
        self.df.dropna(axis=1, how="all", inplace=True)
        self.mapping_dict = {
            "Title": "spectrum_title",
            "File": "raw_data_location",
            "Scan": "spectrum_id",
            "RT(min)": "retention_time_seconds",
            "PrecursorMZ": "exp_mz",
            "Charge": "charge",
            "Peptide": "sequence",
            "Modification": "modifications",
            "PeptideMass": "glyco_decipher:PeptideMass",
            "Protein": "protein_id",
            "GlycoSite": "glyco_decipher:GlycoSite",
            "PeptideScore": "glyco_decipher:PeptideScore",
            "PeptideFDR": "glyco_decipher:PeptideFDR",
            "CoreMatch": "glyco_decipher:CoreMatch",
            "CoreScore": "glyco_decipher:CoreScore",
            "CoreFucosed": "glyco_decipher:CoreFucosed",
            "DiagnosticIon(CoreFucosed)": "glyco_decipher:DiagnosticIon(CoreFucosed)",
            "Bisecting": "glyco_decipher:Bisecting",
            "DiagnosticIon(Bisecting)": "glyco_decipher:DiagnosticIon(Bisecting)",
            "GlycanExpMass": "glyco_decipher:GlycanExpMass",
            "GlycanID": "glyco_decipher:GlycanID",
            "GlycanComposition": "glycan_composition",
            "GlycanMass": "glyco_decipher:GlycanMass",
            "GlycanScore": "glyco_decipher:GlycanScore",
            "GlycanFDR": "glyco_decipher:GlycanFDR",
            "Delta(Da)": "glyco_decipher:Delta(Da)",
            "Delta(ppm)": "glyco_decipher:Delta(ppm)",
        }

        self.df.rename(columns=self.mapping_dict, inplace=True)
        self.df.columns = self.df.columns.str.lstrip(" ")
        # if not "modifications" in self.df.columns:
        # self.df["modifications"] = ""
        self.reference_dict.update({k: None for k in self.mapping_dict.values()})
        self.metadata = {
            "File Origin": "GlycoDecipher",
            "Version": [1.0],
            "bigger_scores_better": True,
            "validation_score_field": "glyco_decipher:PeptideScore",
            "Parser": "pyiohat/parsers/ident/glyco_decipher_1_parser.py",
        }

    @classmethod
    def check_parser_compatibility(cls, file):
        """Assert compatibility between file and parser.

        Args:
            file (str): path to input file

        Returns:
            bool: True if parser and file are compatible

        """
        is_tsv = file.as_posix().endswith(".txt")
        with open(file.as_posix()) as f:
            try:
                head = "".join([next(f) for _ in range(1)])
            except StopIteration:
                head = ""
        head = set(head.rstrip("\n").split("\t"))
        ref_columns = {
            "File",
            "Title",
            "Scan",
            "RT(min)",
            "PrecursorMZ",
            "Charge",
            "Peptide",
            "Modification",
            "PeptideMass",
            "Protein",
            "GlycoSite",
            "PeptideScore",
            "PeptideFDR",
            "CoreMatch",
            "CoreScore",
            "CoreFucosed",
            "DiagnosticIon(CoreFucosed)",
            "Bisecting",
            "DiagnosticIon(Bisecting)",
            "GlycanExpMass",
            "GlycanID",
            "GlycanComposition",
            "GlycanMass",
            "GlycanScore",
            "GlycanFDR",
            "Delta(Da)",
            "Delta(ppm)",
        }
        columns_match = len(ref_columns.difference(head)) == 0
        return is_tsv and columns_match

    def convert_glycan_composition(self):
        """
        Converts the oGlyco formatted glycan compositions into unified glycan
        compositions
        """
        return_df = self.df["glycan_composition"].apply(self.transform_glycan_entry)
        return return_df

    def transform_glycan_entry(self, entry):
        # Find all patterns of <single_letter>(<integer>)
        glyco_lookup = {
            "Hex": "Hex",
            "HexNAc": "HexNAc",
            "Fuc": "dHex",
            "NeuAc": "NulNAcA",
            "NeuGc": "NulNGcA",
        }
        matches = re.findall(r"([A-Za-z0-9]+)\((\d+)\)", entry)

        if not matches:
            return None  # Handle cases with no valid pattern

        transformed_parts = []
        for org, number in matches:
            if org in glyco_lookup:
                # Replace the letter with full name defined in lookup
                transformed_parts.append(f"{glyco_lookup[org]}({number})")
            else:
                raise ValueError(
                    f"Unknown single monosaccharide letter in GlycoDecipher results: '{org}' found in '{entry}'"
                )

        return "".join(transformed_parts)

    def adjust_modifications(self):
        """
        Adjusts the pGlyco modifications (<pos,<name>[aa];) to pyiohat style (<name>:<pos>).
        """

        # Apply the transformation function to the 'modification' column
        return self.df["modifications"].apply(self.transform_mod_entry)

    def transform_mod_entry(self, entry):
        if pd.isna(entry) or not entry:
            return ""

        # Split the string by ';' and filter out any empty strings resulting from
        # a trailing semicolon.
        modifications = [mod for mod in entry.split(";") if mod]

        transformed_mods = []
        for mod in modifications:
            # Use regex to find the position and name
            match = re.search(r"(\d+),(.+)\(.+\)", mod)
            if match:
                pos = match.group(1)
                name = match.group(2)
                transformed_mods.append(f"{name}:{pos}")

        return ";".join(transformed_mods)

    def unify(self):
        """
        Primary method to read and unify engine output.

        Returns:
            self.df (pd.DataFrame): unified dataframe
        """
        self.df["search_engine"] = "glyco_decipher_1"
        self.df["glycan_is_decoy"] = False
        decoy_tag = self.params.get("decoy_tag", "decoy_")
        self.df["peptide_is_decoy"] = self.df["protein_id"].str.contains(decoy_tag)
        self.df["retention_time_seconds"] *= 60.0
        self.df["glycan_composition"] = self.convert_glycan_composition()
        self.df["modifications"] = self.adjust_modifications()
        self.process_unify_style()

        return self.df
