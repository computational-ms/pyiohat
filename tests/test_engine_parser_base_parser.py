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
    flat_precursor_mzs = [
        mz for specs in rt_lookup.values() for mz in specs["precursor_mz"]
    ]
    # 2. Use np.isnan() since np.nan == np.nan evaluates to False in Python
    assert all(np.isnan(mz) for mz in flat_precursor_mzs)
    # check consistency
    assert 2450 in rt_lookup
    assert rt_lookup[2450]["rt"][0] == pytest.approx(1534.4619140625)
    assert np.isnan(rt_lookup[2450]["precursor_mz"])
    rt_index = rt_lookup[2450]["rt"].index(1534.4619140625)
    assert rt_lookup[2450]["lineage_root"][rt_index] == "path/for/glory.mzML"
    # 1. Find the position (index) of this specific retention time
    rt_index = rt_lookup[2450]["rt"].index(1534.4619140625)
    # 2. Verify that the precursor m/z at the exact same position is NaN
    assert rt_lookup[2450]["precursor_mz"][rt_index] is np.nan
