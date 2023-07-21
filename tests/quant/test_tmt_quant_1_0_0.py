#!/usr/bin/env python
import pytest

from pyiohat.parsers.quant.tmt_quant_parser_1_0_0 import (
    TMTQuantParser,
)


def test_engine_parsers_tmt_quant_init():
    input_file = pytest._test_path / "data" / "mapping_data" / "trunc_tmt_quant.csv"
    rt_lookup_path = pytest._test_path / "data" / "TMT16_ursgal_lookup.csv"
    parser = TMTQuantParser(
        input_file,
        params={
            "rt_pickle_name": rt_lookup_path,
        },
    )


def test_engine_parsers_tmt_quant_check_parser_compatibility():
    input_file = pytest._test_path / "data" / "mapping_data" / "trunc_tmt_quant.csv"

    assert TMTQuantParser.check_parser_compatibility(input_file) is True


def test_engine_parsers_tmt_quant_file_not_matches_parser():
    input_file = (
        pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_msfragger_3.tsv"
    )

    assert TMTQuantParser.check_parser_compatibility(input_file) is False


def test_engine_parsers_tmt_quant_unify_row():
    input_file = pytest._test_path / "data" / "mapping_data" / "trunc_tmt_quant.csv"
    rt_lookup_path = pytest._test_path / "data" / "TMT16_ursgal_lookup.csv"

    parser = TMTQuantParser(
        input_file,
        params={
            "rt_pickle_name": rt_lookup_path,
        },
    )
    df = parser.unify()
    assert len(df) == 7
    assert pytest.approx(df["retention_time_seconds"].mean()) == 781.94380
    assert pytest.approx(df["quant_value"].mean()) == 249596.16588
    assert pytest.approx(df["accuracy_mz"].mean()) == 0.009146615862846
