"""Parser handler."""

import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd
from chemical_composition import ChemicalComposition
from loguru import logger
from unimod_mapper.unimod_mapper import UnimodMapper


class BaseParser:
    """Base class of all parser types."""

    def __init__(self, input_file, params, immutable_peptides=None):
        """Initialize parser.

        Args:
            input_file (str): path to input file
            params (dict): ursgal param dict
            immutable_peptides (list, optional): list of immutable peptides
        """
        self.input_file = input_file
        self.immutable_peptides = immutable_peptides
        if params is None:
            params = {}
        self.params = params
        self.xml_file_list = self.params.get("xml_file_list", None)
        self.mod_mapper = UnimodMapper(xml_file_list=self.xml_file_list)
        self.params["mapped_mods"] = self.mod_mapper.map_mods(
            mod_list=self.params.get("modifications", [])
        )
        self.mapped_mod_names = set(
            mod["name"] for mod in self.params["mapped_mods"]["fix"]
        ) | set(mod["name"] for mod in self.params["mapped_mods"]["opt"])
        self.mod_mapper.read_mapped_mods_as_df(self.params["mapped_mods"])
        self.cc = ChemicalComposition(
            unimod_file_list=self.params.get("xml_file_list", None)
        )
        self.style = None

    @classmethod
    def check_parser_compatibility(cls, file):
        """Assert compatibility between file and parser.

        Args:
            file (str): path to input file

        Returns:
            bool: True if parser and file are compatible

        """
        return False

    def _load_model(self, file):
        with open(
            Path(__file__).parent / "models" / "base_parser_model.json", "r"
        ) as jf:
            model_definition = json.load(jf)
        with open(Path(__file__).parent / "models" / file, "r") as jf:
            model_definition.extend(json.load(jf))
        model_definition = {item["name"]: item["dtype"] for item in model_definition}
        return model_definition

    def _read_meta_info_lookup_file(self):
        """Read meta info lookup file.

        Returns:
            rt_lookup (dict of int: dict): dict with spectrum ids as top level key
                                           values are new dicts with all rt values as keys
                                           values are lists with [file, precursor_mz]
        """
        rt_lookup = {}
        with open(self.params["rt_pickle_name"], mode="r") as meta_csv:
            meta_reader = csv.DictReader(meta_csv)
            for row in meta_reader:
                ms_level = row.get("ms_level", "")
                rt = float(row["rt"])
                if row["rt_unit"] == "minute" or row["rt_unit"] == "min":
                    rt *= 60.0
                if row["precursor_mz"] == "":
                    precursor_mz = np.nan
                else:
                    precursor_mz = float(row["precursor_mz"])
                spec_id = int(row["spectrum_id"])
                if spec_id not in rt_lookup:
                    rt_lookup[spec_id] = {
                        "lineage_root": [row["lineage_root"]],
                        "precursor_mz": [precursor_mz],
                        "rt": [rt],
                        "ms_level": [ms_level],
                    }
                else:
                    rt_lookup[spec_id]["lineage_root"].append(row["lineage_root"])
                    rt_lookup[spec_id]["precursor_mz"].append(precursor_mz)
                    rt_lookup[spec_id]["rt"].append(rt)
                    rt_lookup[spec_id]["ms_level"].append(ms_level)
        return rt_lookup

    def sanitize(self):
        """Perform
        .dataframe sanitation steps.

        - Cast defined types
        - Columns that were not filled in but should exist in the unified format are added and set to None
        - Columns in the dataframe which could not be properly mapped are removed (warning is raised)
        Operations are performed inplace on self.df
        """
        # Set missing columns to None and reorder columns in standardized manner
        if self.style and self.style in (
            "pglyco_db_style_1",
            "glyco_decipher_style_1",
            "ptmshepherd_style_1",
        ):
            self.required_headers = self._load_model("ident_glyco_parser_model.json")
        col_order = pd.Series(self.required_headers.keys())
        new_cols = col_order[~col_order.isin(self.df.columns)].to_list()
        self.df.loc[:, new_cols] = pd.NA
        self.df = self.df.loc[
            :,
            col_order.tolist()
            + sorted(self.df.columns[~self.df.columns.isin(col_order)].tolist()),
        ]
        str_cols = [c for c, dtype in self.required_headers.items() if dtype == "str"]
        self.df.loc[:, str_cols] = self.df[str_cols].fillna("")
        self.df = self.df.astype(self.required_headers)
        # Ensure there are not any column that should not be
        if hasattr(self, "mapping_dict"):
            new_cols = set(self.mapping_dict.keys())
        else:
            new_cols = set()
        additional_cols = set(self.df.columns).difference(
            set(self.required_headers.keys()) | new_cols
        )
        unmapped_add_cols = [c for c in additional_cols if ":" not in c]
        if len(unmapped_add_cols) > 0:
            logger.warning(
                f"Some engine level columns ({unmapped_add_cols}) were not properly mapped and removed."
            )
            self.df.drop(columns=unmapped_add_cols, inplace=True, errors="ignore")

        # Drop unwanted duplicated rows
        init_len = len(self.df)
        self.df.drop_duplicates(inplace=True)
        rows_dropped = init_len - len(self.df)
        if rows_dropped != 0:
            logger.warning(
                f"{rows_dropped} duplicated rows were dropped in output csv."
            )
