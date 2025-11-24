"""Engine parser."""

import itertools

import pandas as pd
import regex as re
from loguru import logger

from pyiohat.parsers.ident_base_parser import IdentBaseParser
from pprint import pprint
from itertools import combinations
from chemical_composition import chemical_composition_kb


class Casanovo_5_Parser(IdentBaseParser):
    """File parser for Casanovo 5"""

    def __init__(self, *args, **kwargs):
        """Initialize parser.

        Reads in data file and provides mappings.
        """
        super().__init__(*args, **kwargs)
        self.style = "casanovo_style_5"
        # 15N handling missing for now
        if self.params.get("label", "") == "15N":
            raise NotImplementedError

        self.df = pd.read_csv(self.input_file, delimiter="\t")
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
            "pre": "pre",
            "post": "post",
            "start": "start",
            "end": "end",
            "opt_ms_run[1]_aa_scores": "casanovo:opt_ms_run[1]_aa_scores",
        }
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
            "File Origin": "Casanovo",
            "Version": [4.0, 4.1, 4.2, 4.3],
            "bigger_scores_better": True,
            "validation_score_field": "msfragger:hyperscore",
            "Parser": "pyiohat/parsers/ident/casanovo_5_parser.py",
        }

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
        }
        columns_match = len(ref_columns.difference(head)) == 0
        return is_tsv and columns_match

    def _map_mod_translation(self, row, map_dict):
        """Replace single mod string.

        Args:
            row (str): unprocessed modification string
            map_dict (dict): mod mapping dict

        Returns:
            mod_str (str): formatted modification string
        """
        mod_str = ""
        if row == "" or row == [""]:
            return mod_str
        for mod in row:
            mass = match.group(1) if (match := re.search(r"\(([^)]+)\)", mod)) else None
            if mass == None:
                continue
            # pprint(f"Searching for {mass} in {map_dict}")
            name = map_dict[mass]
            if len(name) > 0:
                for m in name:
                    # if there are digits in mod, process them normally, otherwise check for N-term
                    str_regex_on_mod = re.search(r"^\d+", mod)
                    # pos is 1 so rest of the code flows as expected
                    pos = None
                    # if digits are found, extract them and override position
                    if str_regex_on_mod is not None:
                        pos = int(re.search(r"^\d+", mod).group(0))
                    # TO DO: Does this work if same mod at pos 0 and 1? E.g. TMT
                    if any(
                        [
                            "N-term" in p
                            for p in self.mod_mapper.query(f"`Name` == '{m}'")[
                                "position"
                            ].to_list()
                        ]
                    ) and ((pos == None and "N-term" in mod) or pos == 1):
                        pos = 0
                    else:
                        if pos == None:
                            continue
                        pos = int(re.search(r"^\d+", mod).group(0))
                    mod_str += f"{m}:{pos};"
            else:
                return "NON_MAPPABLE"
        return mod_str

    def translate_mods(self):
        """
        Replace internal modification nomenclature with formatted modification strings.

        Returns:
            (pd.Series): column with formatted mod strings
        """
        self.df["modifications"] = self.df["modifications"].astype(str)
        # print(self.df["modifications"])
        mod_split_col = self.df["modifications"].fillna("").str.split(", ")
        unique_mods = set().union(*mod_split_col.apply(set)).difference({""})
        unique_mod_masses = {
            match.group(1)
            for m in unique_mods
            if (match := re.search(r"\(([^)]+)\)", m))
        }
        # Map single mods
        potential_names = {
            m: [name for name in self.mod_mapper.mass_to_names(float(m), decimals=4)]
            for m in unique_mod_masses
        }
        # Map multiple mods
        for n in [2, 3]:
            for unmapped_mass in {k: v for k, v in potential_names.items() if v == []}:
                potential_mods = [
                    name[1]
                    for name in self.mod_mapper.mass_to_combos(
                        float(unmapped_mass), n=n, decimals=4
                    )
                ]
                if len(potential_mods) == 1:
                    potential_names[unmapped_mass] = potential_mods[0]
        non_mappable_mods = {
            k: len(
                [
                    m
                    for m in list(
                        itertools.chain.from_iterable(
                            mod_split_col.apply(list).to_list()
                        )
                    )
                    if k in m
                ]
            )
            for k, v in potential_names.items()
            if v == []
        }
        non_mappable_percent = pd.Series(
            [v / len(self.df) for v in non_mappable_mods.values()], dtype="float64"
        )
        if any(non_mappable_percent > 0.001):
            raise ValueError(
                "Some modifications found in more than 0.1% of PSMs cannot be mapped."
            )
        if len(non_mappable_percent) > 0:
            logger.warning(
                "Some modifications found in less than 0.1% of PSMs cannot be mapped and were removed."
            )
        # pprint(f"Potential names for modifications reported: {potential_names}")
        mods_translated = mod_split_col.apply(
            self._map_mod_translation, map_dict=potential_names
        )

        return mods_translated.str.rstrip(";")

    def unify(self):
        """
        Primary method to read and unify engine output.

        Returns:
            self.df (pd.DataFrame): unified dataframe
        """
        self.df["search_engine"] = "casanovo_5_0"
        self.df["retention_time_seconds"] *= 60.0
        self.df["exp_mz"] = self._calc_mz(
            mass=self.df["msfragger:precursor_neutral_mass_da"],
            charge=self.df["charge"],
        )
        self.df["modifications"] = self.translate_mods()
        self.df = self.df.loc[
            ~self.df["modifications"].str.contains("NON_MAPPABLE", regex=False), :
        ]
        self.process_unify_style()

        return self.df
