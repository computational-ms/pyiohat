import numpy as np
import pytest

from pyiohat.parsers.base_parser import BaseParser


def test_uninitialized_parser_compatiblity_is_false():
    input_file = (
        pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_xtandem_alanine.xml"
    )
    compat = BaseParser.check_parser_compatibility(input_file)
    assert compat is False


def test_base_parser_read_rt_lookup_file_wo_precursor_mz_info():
    rt_lookup_path = (
        pytest._test_path / "data" / "BSA1_ursgal_lookup_no_precursor_mz.csv"
    )
    input_file = (
        pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_xtandem_alanine.xml"
    )

    bp = BaseParser(input_file, params={"rt_pickle_name": rt_lookup_path})
    rt_lookup = bp._read_meta_info_lookup_file()
    assert len(rt_lookup) == 9
    precursor_mzs = [list(specs.values())[-1][-1] for specs in rt_lookup.values()]
    assert set(precursor_mzs) == {np.nan}
    # check consistency
    assert 2450 in rt_lookup
    assert 1534.4619140625 in rt_lookup[2450]
    assert rt_lookup[2450][1534.4619140625][0] == "path/for/glory.mzML"
    assert rt_lookup[2450][1534.4619140625][1] is np.nan
