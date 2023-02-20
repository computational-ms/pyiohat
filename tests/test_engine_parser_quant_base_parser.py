#!/usr/bin/env python
import pytest

from pyiohat.parsers.quant_base_parser import QuantBaseParser


def test_engine_parsers_QuantBaseParser_init():
    input_file = (
        pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_xtandem_alanine.xml"
    )
    rt_lookup_path = pytest._test_path / "data" / "_ursgal_lookup.csv"
    db_path = pytest._test_path / "data" / "test_Creinhardtii_target_decoy.fasta"

    parser = QuantBaseParser(
        input_file,
        params={
            "cpus": 2,
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "modifications": [
                {
                    "aa": "M",
                    "type": "opt",
                    "position": "any",
                    "name": "Oxidation",
                },
                {
                    "aa": "C",
                    "type": "fix",
                    "position": "any",
                    "name": "Carbamidomethyl",
                },
                {
                    "aa": "*",
                    "type": "opt",
                    "position": "Prot-N-term",
                    "name": "Acetyl",
                },
            ],
        },
    )


def test_engine_parsers_QuantBaseParser_check_parser_compatibility_non_existing():
    # should always return False
    assert QuantBaseParser.check_parser_compatibility("whatever") is False


def test_engine_parsers_QuantBaseParser_check_parser_compatibility_existing():
    # should always return False
    assert (
        QuantBaseParser.check_parser_compatibility(
            pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_xtandem_alanine.xml"
        )
        is False
    )
