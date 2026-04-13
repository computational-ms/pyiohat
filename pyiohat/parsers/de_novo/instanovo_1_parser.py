"""Engine parser."""

import itertools

import pandas as pd
import regex as re
from loguru import logger

from pyiohat.parsers.de_novo_base_parser import DeNovoBaseParser
from pprint import pprint
from itertools import combinations
from chemical_composition import chemical_composition_kb
from pathlib import Path
from unimod_mapper.unimod_mapper import UnimodMapper

class Instanovo_1_Parser(DeNovoBaseParser):
    """File parser for Instanovo 1"""

    def __init__(self, *args, **kwargs):
        """Initialize parser.

        Reads in data file and provides mappings.
        """
        super().__init__(*args, **kwargs)
        self.style = "instanovo_style_1"

        self.mod_mapper = UnimodMapper()

        self.df = pd.read_csv(self.input_file, sep=",")

        self.df.dropna(axis=1, how="all", inplace=True)

        # pprint(f"direct file read from msf out")
        # pprint(self.df)
        self.mapping_dict = {
            "scan_number": "spectrum_id",
            "precursor_mz" : "exp_mz",
            "precursor_charge": "charge",
            "experiment_name":"raw_data_location",
            "spectrum_id": "spectrum_title",
            "retention_time_seconds":"retention_time_seconds",
            "diffusion_predictions_tokenised" : "instanovo:diffusion_predictions_tokenised",
            "diffusion_predictions" : "instanovo:diffusion_predictions",
            "diffusion_log_probabilities" : "instanovo:diffusion_log_probabilities",
            "transformer_predictions":"instanovo:transformer_predictions",
            "transformer_predictions_tokenised":"instanovo:transformer_predictions_tokenised",
            "transformer_log_probabilities":"instanovo:transformer_log_probabilities",
            "transformer_token_log_probabilities":"instanovo:transformer_token_log_probabilities",
            "final_prediction":"sequence",
            "modifications":"modifications",
            "final_prediction_tokenised":"instanovo:final_prediction_tokenised",
            "final_log_probabilities":"instanovo:final_log_probabilities",
            "selected_model":"instanovo:selected_model",
            "precursor_mass_match":"instanovo:precursor_mass_match",
        }
        # pprint(f"mapping dict")
        # pprint(self.mapping_dict)
        self.df.rename(columns=self.mapping_dict, inplace=True)

        # pprint(f"renamed df")
        # pprint(self.df)
        self.df.columns = self.df.columns.str.lstrip(" ")
        self.reference_dict.update({k: None for k in self.mapping_dict.values()})
        self.metadata = {
            "File Origin": "Instanovo",
            "Version": [4.0, 4.1, 4.2, 4.3],
            "bigger_scores_better": True,
            "validation_score_field": "instanovo:transformer_log_probabilities",
            "Parser": "pyiohat/parsers/ident/instanovo_1_parser.py",
        }

    def parse_sequence_modifications(self):
        """
        Parse sequence column to extract mass difference and modifications.
        """
        seq_col = "sequence"

        if seq_col not in self.df.columns:
            raise KeyError(
                f"Column '{seq_col}' not found in DataFrame. Available columns: {list(self.df.columns)}"
            )

        seq_col_index = self.df.columns.get_loc(seq_col)

        # Initialize modifications and mass difference column if it doesn't exist
        if "modifications" not in self.df.columns:
            self.df["modifications"] = None
   
        # Initialize lists for new data
        modifications_list = []
        cleaned_sequences = []

        # Process ALL rows (starting from 0, not 1)
        for idx in range(len(self.df)):
            seq = self.df.iloc[idx, seq_col_index]

            if pd.isna(seq):
                modifications_list.append(None)
                cleaned_sequences.append(seq)
                continue

            seq = str(seq)

            # Extract all modifications like [UNIMOD:35] 
            mod_pattern = r"\[([A-Za-z]+[A-Za-z0-9:.\-()]+)\]"
            modifications = []

            for match in re.finditer(mod_pattern, seq):
                mod_name = match.group(1)
                position = len(
                    re.sub(
                        r"\[([A-Za-z]+[A-Za-z0-9:.\-()]+)\]", "", seq[: match.start()]
                    )
                )
                

                # Look up unimod modifications by ID
                # Try to resolve mod_name as a unimod ID to get canonical name(s)
                try:
                    
                    unimod_id_str = mod_name.split(":")[-1]
                    matched_names = self.mod_mapper.id_to_name(unimod_id_str)
                    if matched_names:
                        modifications.append(f"{matched_names[0]}:{position}")
                    else:
                        # Fallback if the mapper returned an empty list
                        modifications.append(f"{mod_name}:{position}")
                except Exception as e:
                    logger.info(f"DEBUG: Exception during mod lookup for '{mod_name}': {e}")
 

            if modifications:
                modifications_list.append(";".join(modifications))
            else:
                modifications_list.append(None)

            # Remove all modification brackets from sequence
            cleaned_seq = re.sub(r"\[[^\]]+\]", "", seq)
            cleaned_seq = cleaned_seq.lstrip("-")
            cleaned_sequences.append(cleaned_seq)

        # Add Mass Difference column (insert after sequence column)

        # Update modifications and sequence columns
        self.df["modifications"] = modifications_list
        self.df[seq_col] = cleaned_sequences

  

    @classmethod
    def check_parser_compatibility(cls, file):
        """Assert compatibility between file and parser.

        Args:
            file (str): path to input file

        Returns:
            bool: True if parser and file are compatible

        """
        is_instanovo_csv = file.as_posix().endswith(".instanovo.csv")

        with open(file.as_posix()) as f:
            try:
                head = next(f)
            except StopIteration:
                head = ""

        head_raw = head.rstrip("\n")
        delimiter = "\t" if "\t" in head_raw else ","
        head = set(head_raw.split(delimiter))
        ref_columns = {
            "scan_number",
            "precursor_mz",
            "precursor_charge",
            "experiment_name",
            "spectrum_id",
            "diffusion_predictions_tokenised",
            "diffusion_predictions",
            "diffusion_log_probabilities",
            "transformer_predictions",
            "transformer_predictions_tokenised",
            "transformer_log_probabilities",
            "transformer_token_log_probabilities",
            "final_prediction",
            "final_prediction_tokenised",
            "final_log_probabilities",
            "selected_model",
            "precursor_mass_match",
        }
        columns_match = len(ref_columns.difference(head)) == 0
        return is_instanovo_csv and columns_match


    def unify(self):
        """
        Primary method to read and unify engine output.

        Returns:
            self.df (pd.DataFrame): unified dataframe
        """
        self.df["search_engine"] = "instanovo_1_2_2"

    
        self.parse_sequence_modifications()

        # Drop rows where sequence is null (no prediction made)
        self.df = self.df.dropna(subset=["sequence"]).reset_index(drop=True)

        self.process_unify_style()

        return self.df
