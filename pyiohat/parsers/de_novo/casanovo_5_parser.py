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


class Casanovo_5_Parser(DeNovoBaseParser):
    """File parser for Casanovo 5"""

    def __init__(self, *args, **kwargs):
        """Initialize parser.

        Reads in data file and provides mappings.
        """
        super().__init__(*args, **kwargs)
        self.style = "casanovo_style_1"

        self.df = self._read_and_clean_mztab(self.input_file)
        self.df.dropna(axis=1, how="all", inplace=True)

        # pprint(f"direct file read from msf out")
        # pprint(self.df)
        self.mapping_dict = {
            "sequence": "sequence",
            "PSM_ID": "casanovo:PSM_ID",
            "accession": "casanovo:accession",
            "unique": "casanovo:unique",
            "database": "casanovo:database",
            "database_version": "casanovo:database_version",
            "search_engine": "casanovo:search_engine",
            "search_engine_score[1]": "casanovo:search_engine_score[1]",
            "modifications": "modifications",
            "retention_time": "retention_time_seconds",
            "charge": "charge",
            "exp_mass_to_charge": "exp_mz",
            "calc_mass_to_charge": "calc_mz",
            "spectra_ref": "spectrum_id",
            "pre": "sequence_pre_aa",
            "post": "sequence_post_aa",
            "start": "sequence_start",
            "end": "sequence_end",
            "opt_ms_run[1]_aa_scores": "casanovo:opt_ms_run[1]_aa_scores",
            "opt_ms_run[1]_proforma": "casanovo:opt_ms_run[1]_proforma",
        }
        # pprint(f"mapping dict")
        # pprint(self.mapping_dict)
        self.df.rename(columns=self.mapping_dict, inplace=True)

        if "retention_time_seconds" not in self.df.columns:
            self.df["retention_time_seconds"] = None
        # pprint(f"renamed df")
        # pprint(self.df)
        self.df.columns = self.df.columns.str.lstrip(" ")
        if not "modifications" in self.df.columns:
            self.df["modifications"] = ""
        self.reference_dict.update({k: None for k in self.mapping_dict.values()})
        self.metadata = {
            "File Origin": "Casanovo",
            "Version": [4.0, 4.1, 4.2, 4.3],
            "bigger_scores_better": True,
            "validation_score_field": "casanovo:search_engine_score[1]",
            "Parser": "pyiohat/parsers/ident/casanovo_5_parser.py",
        }

    def _read_and_clean_mztab(self, filepath):
        """
        Read mzTab file and get rid of rows that are not needed by pyiohat.

        Args:
             filepath (str or Path): path to mzTab file

        Returns:
            pd.DataFrame: cleaned dataframe
        """
        skip = next(
            i for i, line in enumerate(open(filepath)) if line.startswith("PSH\t")
        )
        df = pd.read_csv(filepath, sep="\t", skiprows=skip, header=0)
        df = df.iloc[:, 1:]
        df = df.reset_index(drop=True)
        return df

    def parse_sequence_modifications(self):
        """
        Parse the mzTab modifications column into pyiohat format.
        mzTab format: '3-Carbamidomethyl (C):UNIMOD:4' or multiple separated by '|'
        pyiohat format: 'Carbamidomethyl:3;Oxidation:6'
        """
        if "modifications" not in self.df.columns:
            self.df["modifications"] = None
            return

        parsed = []
        for raw in self.df["modifications"]:
            if pd.isna(raw) or raw == "null" or raw == "":
                parsed.append(None)
                continue

            mod_strings = []
            # mzTab can have multiple mods separated by '|'
            for entry in str(raw).split(";"):
                entry = entry.strip()
                # Format: '{position}-{mod_name} ({aa}):UNIMOD:{id}'
                # e.g. '3-Carbamidomethyl (C):UNIMOD:4'
                match = re.match(r"^(\d+)-([^:]+?)(?:\s*\([^)]+\))?:UNIMOD:\d+$", entry)
                if match:
                    position = match.group(1)
                    mod_name = match.group(2).strip()
                    mod_strings.append(f"{mod_name}:{position}")
                else:
                    mass_match = re.match(r"^(\d+)-\[([+-]?\d+\.?\d*)\]$", entry)
                    if mass_match:
                        position = mass_match.group(1)
                        mass = mass_match.group(2)
                        try:
                            matched_names = self.mod_mapper.mass_to_names(
                                float(mass), decimals=4
                            )
                            if len(matched_names) > 0:
                                print(
                                    f"Warning: mass {mass} matched {len(matched_names)} "
                                    f"name(s) {matched_names}; only using first match "
                                    f"'{matched_names[0]}'"
                                )
                                mod_strings.append(
                                    f"{matched_names[0]}:{position}"
                                )  # position 0 as placeholder

                        except Exception as e:
                            print(
                                f"Warning: mod_mapper lookup failed for mass {mass}: no UNIMOD match found"
                            )
                    else:
                        # Fallback: keep raw entry so nothing is silently lost
                        mod_strings.append(entry)

            parsed.append(";".join(mod_strings) if mod_strings else None)

        self.df["modifications"] = parsed

    def _extract_scan_number(self):
        """
        Extract scan number from spectra_ref column.

        Converts mzTab format 'ms_run[1]:controllerType=0 controllerNumber=1 scan=2096'
        to just '2096'
        """
        if "spectrum_id" in self.df.columns:
            # Extract scan number from the format: scan=XXXX
            self.df["spectrum_id"] = self.df["spectrum_id"].str.extract(r"scan=(\d+)")[
                0
            ]

        # def _save_to_csv(self, filepath):
        """
        Save the processed DataFrame to a CSV file.

        Args:
            filepath (str or Path): original mzTab file path

        Returns:
            str: path to the saved CSV file
        """
        # Convert Path object to string if necessary

    #    filepath_str = str(filepath)

    # Create output filename by replacing .mzTab or .tsv with .csv
    #   if filepath_str.endswith(".mzTab"):
    #      output_path = filepath_str.replace(".mzTab", "_processed.csv")
    # else:
    #    output_path = filepath_str + ".csv"

    # Save DataFrame to CSV
    # self.df.to_csv(output_path, index=False)

    # logger.info(f"Processed data saved to: {output_path}")

    # return output_path

    @classmethod
    def check_parser_compatibility(cls, file):
        """Assert compatibility between file and parser.

        Args:
            file (str): path to input file

        Returns:
            bool: True if parser and file are compatible

        """
        is_mztab = file.as_posix().endswith(".mztab")

        with open(file.as_posix()) as f:
            head = next((line for line in f if line.startswith("PSH\t")), "")

        head = set(head.rstrip("\n").split("\t"))
        ref_columns = {
            "PSH",
            "sequence",
            "PSM_ID",
            "accession",
            "unique",
            "database",
            "database_version",
            "search_engine",
            "search_engine_score[1]",
            "modifications",
            "retention_time",
            "charge",
            "exp_mass_to_charge",
            "calc_mass_to_charge",
            "spectra_ref",
            "pre",
            "post",
            "start",
            "end",
            "opt_ms_run[1]_aa_scores",
            "opt_ms_run[1]_proforma",
        }
        columns_match = len(ref_columns.difference(head)) == 0
        return is_mztab and columns_match

    def unify(self):
        """
        Primary method to read and unify engine output.

        Returns:
            self.df (pd.DataFrame): unified dataframe
        """
        self.df["search_engine"] = "casanovo_5_0"

        # Parse sequence modifications before column renaming
        # get the raw data location
        # raw_data_file_name = Path(self.input_file).name.replace(".casanovo.mztab", "") #check naming conventions of pyiohat compared to raw files
        # # Remove the Casanovo WID prefix to get the base file name
        # raw_data_file_name = "_".join(raw_data_file_name.split("_")[4:])
        # self.df["raw_data_location"] = raw_data_file_name

        self.parse_sequence_modifications()

        self._extract_scan_number()

        # Save as CSV file
        # self._save_to_csv(self.input_file)

        # score_col = "casanovo:opt_ms_run[1]_aa_scores" #check if redundant
        # if score_col in self.df.columns:
        # Replace commas with semicolons in the score column data
        # self.df[score_col] = self.df[score_col].str.replace(',', ';', regex=False)

        # Convert retention time from minutes to seconds BEFORE process_unify_style
        if "retention_time_seconds" in self.df.columns:
            self.df["retention_time_seconds"] = self.df[
                "retention_time_seconds"
            ].astype(float)

        self.process_unify_style()

        return self.df
