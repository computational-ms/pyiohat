"""Engine parser."""

import pandas as pd
import regex as re
from loguru import logger

from pyiohat.parsers.ident_base_parser import IdentBaseParser
from pprint import pprint


class PGlyco_3_Parser(IdentBaseParser):
    """File parser for pGlyco 3"""

    def __init__(self, *args, **kwargs):
        """Initialize parser.

        Reads in data file and provides mappings.
        """
        super().__init__(*args, **kwargs)
        self.style = "pglyco_db_style_1"

        self.df = pd.read_csv(self.input_file, delimiter="\t")
        self.df.dropna(axis=1, how="all", inplace=True)
        self.mapping_dict = {
            "GlySpec": "pglyco:GlySpec",
            "PepSpec": "pglyco:PepSpec",
            "RawName": "raw_data_location",
            "Scan": "spectrum_id",
            "RT": "retention_time_seconds",
            "PrecursorMH": "pglyco:PrecursorMH",
            "PrecursorMZ": "exp_mz",
            "Charge": "charge",
            "Rank": "rank",
            "Peptide": "sequence",
            "Mod": "modifications",
            "PeptideMH": "pglyco:PeptideMH",
            "Glycan(H,N,A,F)": "pglyco:Glycan(H,N,A,F)",
            "GlycanComposition": "glycan_composition",
            "PlausibleStruct": "pglyco:PlausibleStruct",
            "GlyID": "pglyco:GlyID",
            "GlyFrag": "pglyco:GlyFrag",
            "GlyMass": "pglyco:GlyMass",
            "GlySite": "pglyco:GlySite",
            "TotalScore": "pglyco:TotalScore",
            "PepScore": "pglyco:PepScore",
            "GlyScore": "pglyco:GlyScore",
            "CoreMatched": "pglyco:CoreMatched",
            "MassDeviation": "pglyco:MassDeviation",
            "PPM": "pglyco:PPM",
            "GlyIonRatio": "pglyco:GlyIonRatio",
            "byIonRatio": "pglyco:byIonRatio",
            "czIonRatio": "pglyco:czIonRatio",
            "GlyDecoy": "glycan_is_decoy",
            "PepDecoy": "peptide_is_decoy",
            "Ion_163.06": "pglyco:Ion_163.06",
            "Ion_366.14": "pglyco:Ion_366.14",
            "Ion_204.09": "pglyco:Ion_204.09",
            "Ion_138.05": "pglyco:Ion_138.05",
            "Ion_292.10": "pglyco:Ion_292.10",
            "Ion_274.09": "pglyco:Ion_274.09",
        }
        # self.mapping_dict = {
        #     v: k
        #     for k, v in self.param_mapper.get_default_params(style=self.style)[
        #         "header_translations"
        #     ]["translated_value"].items()
        # }
        # pprint(f"mapping dict")
        # pprint(self.mapping_dict)
        self.df.rename(columns=self.mapping_dict, inplace=True)
        # pprint(f"renamed df")
        # pprint(self.df)
        self.df.columns = self.df.columns.str.lstrip(" ")
        if not "modifications" in self.df.columns:
            self.df["modifications"] = ""
        self.reference_dict.update({k: None for k in self.mapping_dict.values()})
        self.metadata = {
            "File Origin": "pGlyco",
            "Version": [3.0, 3.1],
            "bigger_scores_better": True,
            "validation_score_field": "pglyco:TotalScore",
            "Parser": "pyiohat/parsers/ident/pglyco_3_parser.py",
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
            "GlySpec",
            "PepSpec",
            "RawName",
            "Scan",
            "RT",
            "PrecursorMH",
            "PrecursorMZ",
            "Charge",
            "Rank",
            "Peptide",
            "Mod",
            "PeptideMH",
            "Glycan(H,N,A,F)",
            "GlycanComposition",
            "PlausibleStruct",
            "GlyID",
            "GlyFrag",
            "GlyMass",
            "GlySite",
            "TotalScore",
            "PepScore",
            "GlyScore",
            "CoreMatched",
            "MassDeviation",
            "PPM",
            "GlyIonRatio",
            "byIonRatio",
            "czIonRatio",
            "GlyDecoy",
            "PepDecoy",
            "Ion_163.06",
            "Ion_366.14",
            "Ion_204.09",
            "Ion_138.05",
            "Ion_292.10",
            "Ion_274.09",
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
        pglyco_glyco_lookup = {
            "H": "Hex",
            "N": "HexNAc",
            "F": "dHex",
            "A": "NulNAcA",
            "G": "NulNGcA",
        }
        matches = re.findall(r"([A-Z])\((\d+)\)", entry)

        if not matches:
            return None  # Handle cases with no valid pattern

        transformed_parts = []
        for letter, number in matches:
            if letter in pglyco_glyco_lookup:
                # Replace the letter with full name defined in lookup
                transformed_parts.append(f"{pglyco_glyco_lookup[letter]}({number})")
            else:
                raise ValueError(
                    f"Unknown single monosaccharide letter in pGlyco results: '{letter}' found in '{entry}'"
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
            match = re.search(r"(\d+),(.+)\[.+\]", mod)
            if match:
                pos = match.group(1)
                name = match.group(2)
                transformed_mods.append(f"{name}:{pos}")

        return ";".join(transformed_mods)

    def convert_is_decoy_columns(self):
        """
        Converts 1/0 integer values in 'glycan_is_decoy' and 'peptide_is_decoy'
        columns to True/False boolean values.
        """
        conversion_map = {1: True, 0: False}

        self.df["glycan_is_decoy"] = self.df["glycan_is_decoy"].map(conversion_map)
        self.df["peptide_is_decoy"] = self.df["peptide_is_decoy"].map(conversion_map)

    def unify(self):
        """
        Primary method to read and unify engine output.

        Returns:
            self.df (pd.DataFrame): unified dataframe
        """
        self.df["search_engine"] = "pglyco_3"
        self.df["spectrum_id"] = self.df["pglyco:PepSpec"].str.split(
            ".",
            expand=True,
        )[1]
        self.df["glycan_composition"] = self.convert_glycan_composition()
        self.convert_is_decoy_columns()
        self.df["modifications"] = self.adjust_modifications()
        self.process_unify_style()

        return self.df
