"""Parser handler."""
import json
from pathlib import Path

import pandas as pd
import uparma
from loguru import logger


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
        self.param_mapper = uparma.UParma()
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

    def sanitize(self):
        """Perform dataframe sanitation steps.

        - Cast defined types
        - Columns that were not filled in but should exist in the unified format are added and set to None
        - Columns in the dataframe which could not be properly mapped are removed (warning is raised)
        Operations are performed inplace on self.df
        """
        # Set missing columns to None and reorder columns in standardized manner
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
