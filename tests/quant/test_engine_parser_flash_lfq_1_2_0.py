#!/usr/bin/env python
import pytest

from pyiohat.parsers.quant.flash_lfq_1_2_0_parser import (
    FlashLFQ_1_2_0_Parser,
)


def test_engine_parsers_flashLFQ_init():
    input_file = pytest._test_path / "data" / "flash_lfq_1_2_0_quantified_peaks.tsv"
    rt_lookup_path = pytest._test_path / "data" / "BSA2_ursgal_lookup.csv"
    parser = FlashLFQ_1_2_0_Parser(
        input_file,
        params={
            "rt_pickle_name": rt_lookup_path,
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


def test_engine_parsers_flashLFQ_check_parser_compatibility():
    input_file = pytest._test_path / "data" / "flash_lfq_1_2_0_quantified_peaks.tsv"

    assert FlashLFQ_1_2_0_Parser.check_parser_compatibility(input_file) is True


def test_engine_parsers_flashLFQ_file_not_matches_parser():
    input_file = (
        pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_msfragger_3.tsv"
    )

    assert FlashLFQ_1_2_0_Parser.check_parser_compatibility(input_file) is False


def test_engine_parsers_flashLFQ_unify_row():
    input_file = pytest._test_path / "data" / "flash_lfq_1_2_0_quantified_peaks.tsv"
    rt_lookup_path = pytest._test_path / "data" / "BSA2_ursgal_lookup.csv"

    parser = FlashLFQ_1_2_0_Parser(
        input_file,
        params={
            "rt_pickle_name": rt_lookup_path,
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
    df = parser.unify()
    assert len(df) == 11
    assert pytest.approx(df["flashlfq:ms2_retention_time"].mean()) == 118178.01055184663
    assert pytest.approx(df["quant_value"].mean()) == 306668.4375


def test_engine_parsers_flashLFQ_extract_mods():
    input_file = pytest._test_path / "data" / "flash_lfq_1_2_0_quantified_peaks.tsv"
    rt_lookup_path = pytest._test_path / "data" / "BSA2_ursgal_lookup.csv"

    parser = FlashLFQ_1_2_0_Parser(
        input_file,
        params={
            "rt_pickle_name": rt_lookup_path,
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
    test_sequence = "ELC[Carbamidomethyl]"
    mods = parser.translate_mods(test_sequence)
    assert mods == "Carbamidomethyl:3"

    test_sequence2 = "ELC[Carbamidomethyl]MMMM[Oxidation]"
    mods = parser.translate_mods(test_sequence2)
    assert mods == "Carbamidomethyl:3;Oxidation:7"

    test_sequence3 = "[Acetyl]ELC[Carbamidomethyl]MMMM[Oxidation]"
    mods = parser.translate_mods(test_sequence3)
    assert mods == "Acetyl:0;Carbamidomethyl:3;Oxidation:7"

    # TODO encode C-terminal mods
    # test_sequence4 = "[Acetyl]ELC[Carbamidomethyl]MMMM[Oxidation][TERMINALMOD]"
    # mods = parser.extract_mods(test_sequence4)
    # assert mods == "Acetyl:0;Carbamidomethyl:3;Oxidation:7;TERMINALMOD:8"
