"""Engine parser."""

import itertools

import pandas as pd
import regex as re
from loguru import logger
from importlib import import_module
from pyiohat.parsers.ident_base_parser import IdentBaseParser
from pyiohat.parsers.base_parser import BaseParser
from pprint import pprint
from itertools import combinations
from chemical_composition import chemical_composition_kb
from pathlib import Path


class PTMShepherd_Parser(IdentBaseParser):
    """File parser for PTMShepherd"""

    def __init__(self, *args, **kwargs):
        """Initialize parser.

        Reads in data file and provides mappings.
        """
        super().__init__(*args, **kwargs)
        self.style = "ptmshepherd_style_1"

        self.df = pd.read_csv(self.input_file, delimiter="\t").drop(
            columns="glycan_composition", errors="ignore"
        )
        self.df.dropna(axis=1, how="all", inplace=True)
        self.mapping_dict = {
            "Spectrum File": "raw_data_location",
            "Charge": "charge",
            "Retention": "retention_time_seconds",
            "Spectrum": "spectrum_title",
            "Peptide": "sequence",
            "Assigned Modifications": "modifications",
            "Protein ID": "protein_id",
            "Observed M/Z": "exp_mz",
            "Observed Mass": "exp_mass",
            "Delta Mass": "mass_delta",
            "Total Glycan Composition": "glycan_composition",
            "Glycan Score": "ptm_shepherd:Glycan Score",
            "Glycan q-value": "ptm_shepherd:Glycan q-value",
        }
        self.df.rename(columns=self.mapping_dict, inplace=True)
        self.df.columns = self.df.columns.str.lstrip(" ")
        self.reference_dict.update({k: None for k in self.mapping_dict.values()})
        self.metadata = self._get_metadata()

    def _get_metadata(self):

        metadata = {
            "Version": ["2.0.5"],
            "Parser": "pyiohat/parsers/ident/ptmshepherd_parser.py",
        }

        if not "search_engine" in self.df.columns:
            # This means the parser is running directly after search + ptmshepherd without pyiohat inbetween
            # search engine must have been msfragger, write bigger_score_better and validation_score_field accordinglly
            metadata.update(
                {
                    "validation_score_field": "msfragger:hyperscore",
                    "bigger_scores_better": True,
                }
            )

        else:
            # This means the parser is running after search + pyiohat + ptmshepherd
            # In this case bigger_score_better and validation_score_field needs to be retrived from the parser of the relavent search engine
            parsers_dict = {
                "xtandem_": "xtandem_alanine",
                "omssa_2_1_9": "omssa_2_1_9_parser",
                "msgfplus_": "msgfplus_2021_03_22_parser",
                "msfragger_4_2": "msfragger_4_parser",
                "msfragger_3_0": "msfragger_3_parser",
                "msamanda_2_0_0_17442": "msamanda_2_parser",
                "mascot_": "mascot_2_6_2_parser",
                "comet_": "comet_2020_01_4_parser",
                "pglyco_3": "pglyco_3_parser",
                "glyco_decipher_1": "glyco_decipher_1_parser",
            }
            search_engine = self.df["search_engine"][1]
            for k, v in parsers_dict.items():
                if k in search_engine:
                    parser_name = v
            original_parser_module = f"{__package__}.{parser_name}"
            import_module(original_parser_module)
            parser_classes = []
            for cat in BaseParser.__subclasses__():
                parser_classes.extend(
                    [
                        c
                        for c in cat.__subclasses__()
                        if c.__module__ == original_parser_module
                    ]
                )
            ParserClass = parser_classes[0]
            parser_instance = ParserClass(
                input_file=self.input_file,
                params=self.params,
                immutable_peptides=self.immutable_peptides,
            )

            original_metadata = parser_instance.metadata

            metadata.update(
                {
                    "validation_score_field": original_metadata[
                        "validation_score_field"
                    ],
                    "bigger_scores_better": original_metadata["bigger_scores_better"],
                    "File Origin": original_metadata["File Origin"],
                }
            )

        return metadata

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
            "Spectrum File",
            "Charge",
            "Retention",
            "Spectrum",
            "Peptide",
            "Assigned Modifications",
            "Protein ID",
            "Observed M/Z",
            "Observed Mass",
            "Calculated M/Z",
            "Calculated Peptide Mass",
            "Delta Mass",
            "MSFragger Localization",
            "Modified Peptide",
            "Observed Modifications",
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

            pos = None
            str_regex_on_mod = re.search(r"^\d+", mod)
            if str_regex_on_mod is not None:
                pos = int(str_regex_on_mod.group(0))

            # Check for N-term in the raw string itself (independent of mapping)
            if pos is None and "N-term" in mod:
                pos = 0
            elif pos is None:
                # If no numeric or N-term position is found, we can't process it.
                continue

            name = map_dict[mass]
            if len(name) > 0:
                for m in name:
                    if any(
                        [
                            "N-term" in p
                            for p in self.mod_mapper.query(f"`Name` == '{m}'")[
                                "position"
                            ].to_list()
                        ]
                    ) and ((pos == None and "N-term" in mod) or pos == 1):
                        pos = 0

                    mod_str += f"{m}:{pos};"
            else:
                mod_str += f"NON_MAPPABLE:{pos};"
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
        # print(f"potential_names right after initilizing: {potential_names}")
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

        for unmapped_mass in {k for k, v in potential_names.items() if v == []}:
            mask = (
                self.df["glycan_composition"]
                .astype(str)
                .str.contains(rf"%\s*{unmapped_mass}\b", regex=True, na=False)
            )
            if mask.any():
                matching_mods = (
                    self.df.loc[mask, "glycan_composition"].str.split(", ").explode()
                )
                matching_prefixes = {
                    mod.split(" % ")[0]
                    for mod in matching_mods
                    if re.search(rf"%\s*{unmapped_mass}\b", mod)
                }
                if len(matching_prefixes) == 1:
                    potential_names[unmapped_mass] = [matching_prefixes.pop()]
                elif len(matching_prefixes) > 1:
                    matching_prefixes_list = list(matching_prefixes)
                    logger.warning(
                        f"Ambiguous mapping for mass {unmapped_mass}: {matching_prefixes}\n{matching_prefixes_list[0]} will be used"
                    )
                    potential_names[unmapped_mass] = [matching_prefixes_list[0]]

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
            logger.warning(
                f"Some modifications found in {non_mappable_percent * 100}% of PSMs cannot be mapped."
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

    def translate_glycans(self):
        """
        Transforms the 'glycan_composition' column based on:
          - Replacing "No Glycan Assigned" with empty string
          - Removing 'Decoy_' prefix if present
          - Stripping content after the first space
          - Mapping each character in the remaining string via monosaccharide_dict

        Returns:
            pandas.Series: processed values aligned to self.df index
        """
        monosaccharide_dict = self.params.get("monosaccharide_dict", {"Fuc": "dHex"})
        s = self.df["glycan_composition"].astype(str)
        # Map the No Glycan Matched case to a blank string
        result = s.replace(["No Glycan Matched", "nan"], "")
        # Remove everything after the first space eg. " % 1702.5814"
        result = result.str.split(" ").str[0]

        def repl(match):
            name = match.group(1)  # e.g. "HexNAc"
            count = match.group(2)  # e.g. "2"
            mapped = monosaccharide_dict.get(name, name)
            return f"{mapped}({count})"

        # Find monomer(amount) pattern and map monomer into generic monomer names
        pattern = re.compile(r"([A-Za-z0-9]+)\((\d+)\)")
        mapped_col = result.apply(lambda raw: pattern.sub(repl, raw))

        return mapped_col

    def peptide_is_decoy(self):
        """
        Returns a boolean Series indicating whether each row should be flagged as a decoy peptide.
        Logic:
          - If 'proteins' exists: check if it contains 'decoy_' anywhere.
          - Else if 'protein_id' exists: check that column instead.
          - Raises KeyError if neither column exists.
        """
        df = self.df
        substring = "decoy_"

        if "proteins" in df.columns:
            # Vectorized substring search; treat NaNs as False
            result = (
                df["proteins"].astype(str).str.contains(substring, case=False, na=False)
            )
        elif "protein_id" in df.columns:
            result = (
                df["protein_id"]
                .astype(str)
                .str.contains(substring, case=False, na=False)
            )
        else:
            raise KeyError(
                "Neither 'proteins' nor 'protein_id' column found in DataFrame."
            )

        return result

    def glycan_is_decoy(self):
        """
        Returns a boolean Series indicating which rows in 'glycan_composition'
        start with the prefix 'Decoy_'.
        It also cleans the DataFrame in place by stripping that prefix where present.
        """
        prefix = "Decoy_"
        col = self.df["glycan_composition"]
        # Calcualte glycan_is_decoy values
        is_decoy = col.str.startswith(prefix)
        # Clean glycan_composition by removing "Decoy_" prefix
        self.df["glycan_composition"] = col.str.removeprefix(prefix)

        return is_decoy

    def unify(self):
        """
        Primary method to read and unify engine output.

        Returns:
            self.df (pd.DataFrame): unified dataframe
        """
        self.df["validation_engine"] = "ptmshepherd"
        self.df["modifications"] = self.translate_mods()
        self.df["glycan_composition"] = self.translate_glycans()
        self.df["peptide_is_decoy"] = self.peptide_is_decoy()
        self.df["glycan_is_decoy"] = self.glycan_is_decoy()

        self.process_unify_style()

        return self.df
