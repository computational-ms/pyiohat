#!/usr/bin/env python
import pandas as pd
import pytest

from pyiohat.parsers.ident.xtandem_alanine import (
    XTandemAlanine_Parser,
)


def test_engine_parsers_xtandem_init():
    input_file = (
        pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_xtandem_alanine.xml"
    )

    parser = XTandemAlanine_Parser(
        input_file,
        params={
            "cpus": 2,
            "enzyme": "(?<=[KR])(?![P])",
            "terminal_cleavage_site_integrity": "any",
            "validation_score_field": {"xtandem_alanine": "x!tandem:hyperscore"},
            "bigger_scores_better": {"xtandem_alanine": True},
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


def test_engine_parsers_xtandem_file_matches_xtandem_parser():
    input_file = (
        pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_xtandem_alanine.xml"
    )

    assert XTandemAlanine_Parser.check_parser_compatibility(input_file) is True


def test_engine_parsers_msfragger_file_not_matches_xtandem_parser():
    input_file = (
        pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_msfragger_3.tsv"
    )

    assert XTandemAlanine_Parser.check_parser_compatibility(input_file) is False


def test_engine_parsers_xtandem_check_dataframe_integrity():
    input_file = (
        pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_xtandem_alanine.xml"
    )
    rt_lookup_path = pytest._test_path / "data" / "_ursgal_lookup.csv"
    db_path = pytest._test_path / "data" / "test_Creinhardtii_target_decoy.fasta"

    parser = XTandemAlanine_Parser(
        input_file,
        params={
            "cpus": 2,
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "enzyme": "(?<=[KR])(?![P])",
            "terminal_cleavage_site_integrity": "any",
            "validation_score_field": {"xtandem_alanine": "x!tandem:hyperscore"},
            "bigger_scores_better": {"xtandem_alanine": True},
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
    assert len(df) == 79
    assert (df["raw_data_location"] == "path/for/glory.mzML").all()
    assert pytest.approx(df["ucalc_mz"].mean()) == 796.4324
    assert pytest.approx(df["exp_mz"].mean()) == 796.71967

    assert df["modifications"].str.contains("Acetyl:0").sum() == 1
    assert df["modifications"].str.contains("Oxidation:").sum() == 23
    assert (
        df["modifications"].str.count("Carbamidomethyl:")
        == df["sequence"].str.count("C")
    ).all()
    assert df["modifications"].str.count(":").sum() == 50


def test_engine_parsers_xtandem_nterminal_mod():
    input_file = (
        pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_xtandem_alanine.xml"
    )
    rt_lookup_path = pytest._test_path / "data" / "_ursgal_lookup.csv"
    db_path = pytest._test_path / "data" / "test_Creinhardtii_target_decoy.fasta"

    parser = XTandemAlanine_Parser(
        input_file,
        params={
            "cpus": 2,
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "enzyme": "(?<=[KR])(?![P])",
            "terminal_cleavage_site_integrity": "any",
            "validation_score_field": {"xtandem_alanine": "x!tandem:hyperscore"},
            "bigger_scores_better": {"xtandem_alanine": True},
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
            "raw_file_location": "test_Creinhardtii_QE_pH11.mzML",
            "15N": False,
        },
    )
    df = parser.unify()
    relevant_row = df[df["sequence"] == "WGLVSSELQTSEAETPGLK"]
    assert relevant_row["modifications"].tolist() == ["Acetyl:0"]


def test_engine_parsers_xtandem_multiple_psms():
    input_file = pytest._test_path / "data" / "multiple_psms_xtandem.xml"
    rt_lookup_path = pytest._test_path / "data" / "_ursgal_lookup.csv"
    db_path = pytest._test_path / "data" / "human_ecoli_test_target_decoy.fasta"

    parser = XTandemAlanine_Parser(
        input_file,
        params={
            "cpus": 2,
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "enzyme": "(?<=[KR])(?![P])",
            "terminal_cleavage_site_integrity": "any",
            "validation_score_field": {"xtandem_alanine": "x!tandem:hyperscore"},
            "bigger_scores_better": {"xtandem_alanine": True},
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
            "raw_file_location": "test_Creinhardtii_QE_pH11.mzML",
            "15N": False,
        },
    )

    # Test file with:
    #    - one sequence in first group
    #    - 3 sequences in second group

    df = parser.unify()
    assert len(df) == 4

    assert set(df["sequence"]) == {
        "ITIPITLRMLIAK",
        "SMMNGGSSPESDVGTDNK",
        "SMMNGGSSPESDVGTDNK",
        "SMMNGGSSPESDVGTDNK",
    }
    assert set(df["spectrum_id"]) == {12833, 14525}
    assert set(df["modifications"]) == {
        "Acetyl:0",
        "",
        "Oxidation:2",
        "Oxidation:3",
    }


def test_engine_parsers_xtandem_map_mod_names():
    input_file = (
        pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_xtandem_alanine.xml"
    )

    parser = XTandemAlanine_Parser(
        input_file,
        params={
            "cpus": 2,
            "enzyme": "(?<=[KR])(?![P])",
            "terminal_cleavage_site_integrity": "any",
            "validation_score_field": {"xtandem_alanine": "x!tandem:hyperscore"},
            "bigger_scores_better": {"xtandem_alanine": True},
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
    test_df = pd.DataFrame({"modifications": [["57.021464:0"]], "sequence": ["CERK"]})
    assert parser.map_mod_names(test_df)["modifications"][0] == "Carbamidomethyl:1"


def test_engine_parsers_xtandem_map_mod_names_nterm():
    input_file = (
        pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_xtandem_alanine.xml"
    )

    parser = XTandemAlanine_Parser(
        input_file,
        params={
            "cpus": 2,
            "enzyme": "(?<=[KR])(?![P])",
            "terminal_cleavage_site_integrity": "any",
            "validation_score_field": {"xtandem_alanine": "x!tandem:hyperscore"},
            "bigger_scores_better": {"xtandem_alanine": True},
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

    row = pd.DataFrame(
        {"modifications": [["57.021464:0", "42.010565:0"]], "sequence": ["CERK"]}
    )
    assert set(parser.map_mod_names(row)["modifications"][0].split(";")) == {
        "Carbamidomethyl:1",
        "Acetyl:0",
    }


def test_engine_parsers_xtandem_get_spec_records():
    input_file = (
        pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_xtandem_alanine.xml"
    )
    parser = XTandemAlanine_Parser(input_file=input_file, params=None)
    spec_records, search_engine = parser.get_spec_records()

    assert len(spec_records) == 79
    assert spec_records[0] == {
        "x!tandem:delta": "0.0057",
        "x!tandem:hyperscore": "14.2",
        "x!tandem:nextscore": "8.0",
        "x!tandem:y_score": "9.9",
        "x!tandem:y_ions": "5",
        "x!tandem:b_score": "0.0",
        "x!tandem:b_ions": "0",
        "sequence": "DDVHNMGADGIR",
        "modifications": ["15.99492:5"],
        "calc_mz": "1315.5700",
        "spectrum_title": "test_Creinhardtii_QE_pH11.10381.10381.3",
        "spectrum_id": "10381",
        "charge": "3",
        "retention_time_seconds": "1943.05878",
    }
    assert spec_records[70] == {
        "x!tandem:delta": "2.022",
        "x!tandem:hyperscore": "10.9",
        "x!tandem:nextscore": "8.0",
        "x!tandem:y_score": "5.0",
        "x!tandem:y_ions": "3",
        "x!tandem:b_score": "7.4",
        "x!tandem:b_ions": "1",
        "sequence": "FHVEALMHADLHFTRKHSR",
        "modifications": ["15.99492:6"],
        "calc_mz": "2348.183",
        "spectrum_title": "test_Creinhardtii_QE_pH11.41649.41649.2",
        "spectrum_id": "41649",
        "charge": "2",
        "retention_time_seconds": "5242.7709",
    }
    assert search_engine == "xtandem_alanine"
