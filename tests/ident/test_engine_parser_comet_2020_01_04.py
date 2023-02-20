#!/usr/bin/env python
import pytest

from pyiohat.parsers.ident.comet_2020_01_4_parser import (
    Comet_2020_01_4_Parser,
)


def test_engine_parsers_comet_init():
    input_file = pytest._test_path / "data" / "BSA1_comet_2020_01_4.mzid"

    parser = Comet_2020_01_4_Parser(
        input_file,
        params={
            "cpus": 2,
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


def test_engine_parsers_comet_check_parser_compatibility():
    msgf_parser_class = Comet_2020_01_4_Parser
    input_file = pytest._test_path / "data" / "BSA1_comet_2020_01_4.mzid"
    assert msgf_parser_class.check_parser_compatibility(input_file) is True


def test_engine_parsers_comet_check_parser_compatibility_fail_with_omssa_file():
    msgf_parser_class = Comet_2020_01_4_Parser
    input_file = (
        pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_omssa_2_1_9.csv"
    )
    assert msgf_parser_class.check_parser_compatibility(input_file) is False


def test_engine_parsers_comet_check_dataframe_integrity():
    input_file = pytest._test_path / "data" / "BSA1_comet_2020_01_4.mzid"
    rt_lookup_path = pytest._test_path / "data" / "BSA1_ursgal_lookup.csv"
    db_path = pytest._test_path / "data" / "BSA.fasta"

    parser = Comet_2020_01_4_Parser(
        input_file,
        params={
            "cpus": 2,
            "enzyme": "(?<=[KR])(?![P])",
            "terminal_cleavage_site_integrity": "any",
            "validation_score_field": {"comet_2020_01_4": "comet:e_value"},
            "bigger_scores_better": {"comet_2020_01_4": False},
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
    df = parser.unify()
    assert pytest.approx(df["ucalc_mz"].mean()) == 457.85944
    assert pytest.approx(df["exp_mz"].mean()) == 457.87625

    assert df["modifications"].str.contains("Acetyl:0").sum() == 5
    assert df["modifications"].str.contains("Oxidation:").sum() == 0
    assert (
        df["modifications"].str.count("Carbamidomethyl:")
        == df["sequence"].str.count("C")
    ).all()
    assert df["modifications"].str.count(":").sum() == 38
    assert (df["raw_data_location"] == "path/for/glory.mzML").all()


def test_engine_parsers_comet_get_version():
    input_file = pytest._test_path / "data" / "BSA1_comet_2020_01_4.mzid"
    parser = Comet_2020_01_4_Parser(input_file, params=None)
    version = parser.get_version()
    assert version == "comet_2020_01_4"


def test_engine_parser_comet_map_mod_mass():
    input_file = pytest._test_path / "data" / "BSA1_comet_2020_01_4.mzid"
    parser = Comet_2020_01_4_Parser(input_file=input_file, params=None)
    fixed_mods, modifications = parser.map_mod_mass()
    assert fixed_mods == {"C": "Carbamidomethyl"}
    assert modifications == {
        "57.021464": "Carbamidomethyl",
        "15.994915": "Oxidation",
        "42.010565": "Acetyl",
    }


def test_engine_parsers_comet_get_peptide_lookup():
    input_file = pytest._test_path / "data" / "BSA1_comet_2020_01_4.mzid"
    parser = Comet_2020_01_4_Parser(input_file=input_file, params=None)
    parser.fixed_mods, parser.mod_mass_map = parser.map_mod_mass()
    peptide_lookup = parser.get_peptide_lookup()
    assert len(peptide_lookup) == 24
    assert "EACFAVEGPK;10:42.010565;" in peptide_lookup.keys()
    assert peptide_lookup["LVTDLTK;"] == ("", "LVTDLTK")


def test_engine_parsers_comet_get_spec_records():
    input_file = pytest._test_path / "data" / "BSA1_comet_2020_01_4.mzid"
    parser = Comet_2020_01_4_Parser(input_file=input_file, params=None)
    parser.fixed_mods, parser.mod_mass_map = parser.map_mod_mass()
    parser.peptide_lookup = parser.get_peptide_lookup()
    spec_records = parser.get_spec_records()

    assert len(spec_records) == 60
    assert spec_records[4] == {
        "comet:num_matched_ions": "3",
        "comet:num_unmatched_ions": "33",
        "comet:xcorr": "0.3634",
        "comet:deltacn": "1.0000",
        "comet:score": "6.4000",
        "comet:e_value": "3.35E+01",
        "charge": "3",
        "sequence": "ECCDKPLLEK",
        "exp_mz": "431.205597",
        "calc_mz": "431.205546",
        "modifications": "Carbamidomethyl:2;Carbamidomethyl:3",
        "spectrum_id": "2569",
    }
    assert spec_records[58] == {
        "comet:num_matched_ions": "2",
        "comet:num_unmatched_ions": "50",
        "comet:xcorr": "0.4654",
        "comet:deltacn": "1.0000",
        "comet:score": "5.7000",
        "comet:e_value": "2.96E+01",
        "charge": "3",
        "sequence": "VPQVSTPTLVEVSR",
        "exp_mz": "504.618805",
        "calc_mz": "504.619112",
        "modifications": "",
        "spectrum_id": "3547",
    }


def test_engine_parsers_comet_map_mods_sequences():
    parser = Comet_2020_01_4_Parser(input_file=None, params=None)
    parser.fixed_mods = {"C": "Carbamidomethyl", "O": "Oxidation"}
    sequence = "ABCDEFCABC"
    mods = parser.map_mods_sequences(sequence)
    assert mods == "Carbamidomethyl:3;Carbamidomethyl:7;Carbamidomethyl:10"

    sequence = "ABCDOEFOCABC"
    mods = parser.map_mods_sequences(sequence)
    assert (
        mods
        == "Carbamidomethyl:3;Carbamidomethyl:9;Carbamidomethyl:12;Oxidation:5;Oxidation:8"
    )

    sequence = "OOOOOCCCCCMMMMM"
    mods = parser.map_mods_sequences(sequence)
    assert (
        mods
        == "Carbamidomethyl:6;Carbamidomethyl:7;Carbamidomethyl:8;Carbamidomethyl:9;Carbamidomethyl:10;Oxidation:1;Oxidation:2;Oxidation:3;Oxidation:4;Oxidation:5"
    )

    sequence = "CD"
    mods = parser.map_mods_sequences(sequence)
    assert mods == "Carbamidomethyl:1"

    sequence = "C"
    mods = parser.map_mods_sequences(sequence)
    assert mods == "Carbamidomethyl:1"

    parser.fixed_mods = {}
    sequence = "ASLDPOCSADK"
    with pytest.raises(ValueError):
        parser.map_mods_sequences(sequence)
