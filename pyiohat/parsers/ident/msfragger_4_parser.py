"""Engine parser."""

import itertools

import pandas as pd
import regex as re
from loguru import logger

from pyiohat.parsers.ident_base_parser import IdentBaseParser
from pprint import pprint
from itertools import combinations
from chemical_composition import chemical_composition_kb


class MSFragger_4_Parser(IdentBaseParser):
    """File parser for MSFragger 4"""

    def __init__(self, *args, **kwargs):
        """Initialize parser.

        Reads in data file and provides mappings.
        """
        super().__init__(*args, **kwargs)
        self.style = "msfragger_style_4"
        # 15N handling missing for now
        if self.params.get("label", "") == "15N":
            raise NotImplementedError

        self.df = pd.read_csv(self.input_file, delimiter="\t")
        self.df.dropna(axis=1, how="all", inplace=True)
        # pprint(f"direct file read from msf out")
        # pprint(self.df)
        self.mapping_dict = {
            v: k
            for k, v in self.param_mapper.get_default_params(style=self.style)[
                "header_translations"
            ]["translated_value"].items()
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
            "File Origin": "MSFragger",
            "Version": [4.0, 4.1, 4.2, 4.3],
            "bigger_scores_better": True,
            "validation_score_field": "msfragger:hyperscore",
            "Parser": "pyiohat/parsers/ident/msfragger_4_parser.py",
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
            "scannum",
            "peptide",
            "charge",
            "peptide_prev_aa",
            "peptide_next_aa",
            "proteins",
            "modification_info",
            "retention_time",
            "precursor_neutral_mass",
            "calc_neutral_pep_mass",
            "hit_rank",
            "massdiff",
            "num_matched_ions",
            "tot_num_ions",
            "hyperscore",
            "nextscore",
            "num_tol_term",
            "num_missed_cleavages",
            "expectscore",
            "best_locs",
            "score_without_delta_mass",
            "best_score_with_delta_mass",
            "second_best_score_with_delta_mass",
            "delta_score",
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

    def annotate_delta_mass(self):

        # glycan_dict = {23.58839: "Hex(2)[23.58839]", 79.24299: "Hex(2)[23.58839];ACDS(12)[55.6546]"}
        # loop through u run dict and collect all glycans and their mass
        glycans = {}
        for mod in self.params["mapped_mods"]["opt"]:
            if "labile" not in mod.keys():
                continue
            glycans[float(mod["mass"])] = f"{mod['name']}[{mod['mass']}]"
        if glycans == {}:
            return [None] * len(self.df["modifications"])
        masses = list(glycans.keys())
        mass_combos = []
        for i in range(1, len(masses) + 1):
            for combo in combinations(masses, i):
                mass_combos.append(list(combo))

        glycan_dict = {}
        for combo in mass_combos:
            value = ""
            key = 0.0
            for m in combo:
                value = value + f"{glycans[m]};"
                key = key + m
            glycan_dict[key] = value
        pprint(glycan_dict)
        annotated_delta_mass = [
            self._map_delta_mass(delta_mass, pep_mass, glycan_dict)
            for delta_mass, pep_mass in zip(
                self.df["mass_difference"], self.df["msfragger:neutral_mass_of_peptide"]
            )
        ]
        return annotated_delta_mass

    def _map_delta_mass(self, delta_mass, pep_mass, mass_glycan_lookup):

        mass_diff = float(delta_mass)
        pep_mass = float(pep_mass)
        if -2 <= round(mass_diff) <= 2:
            return None
        n = 0
        t_mass_diff_L, t_mass_diff_U = self._transform_mass_add_error(
            mass_diff, pep_mass
        )
        for potential_mass in mass_glycan_lookup.keys():
            if t_mass_diff_L <= potential_mass <= t_mass_diff_U:
                return mass_glycan_lookup[potential_mass]
        while True:
            n += 1
            mass_diff = mass_diff - chemical_composition_kb.PROTON
            if n == 4:
                pprint("Give up ----------------------------------------------------")
                pprint(f"delta_mass: {delta_mass}, peptide_mass: {pep_mass}")
                pprint(
                    f"range searched: {self._transform_mass_add_error(mass_diff + 4 * (chemical_composition_kb.PROTON), pep_mass)[0]} to {self._transform_mass_add_error(mass_diff + 4 * (chemical_composition_kb.PROTON), pep_mass)[1]}"
                )
                pprint(
                    f"range searched (1): {self._transform_mass_add_error(mass_diff + 3 * (chemical_composition_kb.PROTON), pep_mass)[0]} to {self._transform_mass_add_error(mass_diff + 3 * (chemical_composition_kb.PROTON), pep_mass)[1]}"
                )
                pprint(
                    f"range searched (2): {self._transform_mass_add_error(mass_diff + 2 * (chemical_composition_kb.PROTON), pep_mass)[0]} to {self._transform_mass_add_error(mass_diff + 2 * (chemical_composition_kb.PROTON), pep_mass)[1]}"
                )
                pprint(
                    f"range searched (3): {self._transform_mass_add_error(mass_diff + (chemical_composition_kb.PROTON), pep_mass)[0]} to {self._transform_mass_add_error(mass_diff + (chemical_composition_kb.PROTON), pep_mass)[1]}"
                )
                pprint(
                    "-----------------------------------------------------------------------"
                )
                return "n=4"
            t_mass_diff_L, t_mass_diff_U = self._transform_mass_add_error(
                mass_diff, pep_mass
            )
            for potential_mass in mass_glycan_lookup.keys():
                if t_mass_diff_L <= potential_mass <= t_mass_diff_U:
                    return mass_glycan_lookup[potential_mass]
        return None

    def _transform_mass_add_error(self, mass, pep_mass):
        if self.params["precursor_mass_tolerance_unit"] == "ppm":
            lower_mass = (
                mass
                - 2
                * self.params["precursor_mass_tolerance_minus"]
                * (mass + pep_mass)
                / 1e6
            )
            upper_mass = (
                mass
                + 2
                * self.params["precursor_mass_tolerance_plus"]
                * (mass + pep_mass)
                / 1e6
            )
        elif self.params["precursor_mass_tolerance_unit"] != "da":
            lower_mass = (mass + pep_mass) - 2 * self.params[
                "precursor_mass_tolerance_minus"
            ]
            upper_mass = (mass + pep_mass) + 2 * self.params[
                "precursor_mass_tolerance_plus"
            ]
        else:
            print(
                "[ERROR] mass tolerance unit {0} not supported".format(
                    self.params["precursor_mass_tolerance_unit"]
                )
            )
            sys.exit(1)
        return lower_mass, upper_mass

    def unify(self):
        """
        Primary method to read and unify engine output.

        Returns:
            self.df (pd.DataFrame): unified dataframe
        """
        self.df["search_engine"] = "msfragger_4_2"
        self.df["retention_time_seconds"] *= 60.0
        self.df["exp_mz"] = self._calc_mz(
            mass=self.df["msfragger:precursor_neutral_mass_da"],
            charge=self.df["charge"],
        )
        self.df["modifications"] = self.translate_mods()
        self.df["annotated_delta_mass"] = self.annotate_delta_mass()
        self.df = self.df.loc[
            ~self.df["modifications"].str.contains("NON_MAPPABLE", regex=False), :
        ]
        self.process_unify_style()

        return self.df
