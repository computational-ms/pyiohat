from pathlib import Path

import pytest

from unify_idents.engine_parsers.ident.msfragger_3_parser import MSFragger_3_Parser


def test_engine_parsers_msfragger_init():
    input_file = (
        pytest._test_path
        / "data"
        / "test_Creinhardtii_QE_pH11_mzml2mgf_0_0_1_msfragger_3.tsv"
    )
    rt_lookup_path = pytest._test_path / "data" / "_ursgal_lookup.csv.bz2"
    db_path = (
        pytest._test_path / "data" / "test_Creinhardtii_target_decoy.fasta"
    )

    parser = MSFragger_3_Parser(
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
            # "omssa_mod_dir": pytest._test_path / "data",
        },
    )


def test_engine_parsers_msfragger_check_parser_compatibility():
    input_file = (
        pytest._test_path
        / "data"
        / "test_Creinhardtii_QE_pH11_mzml2mgf_0_0_1_msfragger_3.tsv"
    )
    assert MSFragger_3_Parser.check_parser_compatibility(input_file) is True


def test_engine_parsers_msfragger_check_dataframe_integrity():
    input_file = (
        pytest._test_path
        / "data"
        / "test_Creinhardtii_QE_pH11_mzml2mgf_0_0_1_msfragger_3.tsv"
    )
    rt_lookup_path = pytest._test_path / "data" / "_ursgal_lookup.csv.bz2"
    db_path = (
        pytest._test_path / "data" / "test_Creinhardtii_target_decoy.fasta"
    )

    parser = MSFragger_3_Parser(
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
            "Raw data location": "/Users/cellzome/Dev/Gits/Ursgal/ursgal_master/example_data/test_Creinhardtii_QE_pH11.mzML",
            "15N": False,
        },
    )
    df = parser.unify()
    assert len(df) == 62
    assert pytest.approx(df["uCalc m/z"].mean()) == 743.15216


def test_engine_parsers_msfragger_merge_mods():
    input_file = pytest._test_path / "data" / "msfragger_merged_mods.tsv"
    rt_lookup_path = pytest._test_path / "data" / "_ursgal_lookup.csv.bz2"
    db_path = (
        pytest._test_path / "data" / "test_Creinhardtii_target_decoy.fasta"
    )

    parser = MSFragger_3_Parser(
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
            "Raw data location": "/Users/cellzome/Dev/Gits/Ursgal/ursgal_master/example_data/test_Creinhardtii_QE_pH11.mzML",
            "15N": False,
        },
    )
    df = parser.unify()
    row = df.iloc[0]
    assert row["Modifications"] == "Acetyl:0;Carbamidomethyl:1"
    assert row["Sequence"] == "CGFSTVGSGFGSR"
    assert (
        row["Raw data location"]
        == "/Users/cellzome/Dev/Gits/Ursgal/ursgal_master/example_data/test_Creinhardtii_QE_pH11.mzML"
    )


def test_engine_parsers_msfragger_single_mods():
    input_file = pytest._test_path / "data" / "msfragger_no_mods.tsv"
    rt_lookup_path = pytest._test_path / "data" / "_ursgal_lookup.csv.bz2"
    db_path = (
        pytest._test_path / "data" / "test_Creinhardtii_target_decoy.fasta"
    )

    parser = MSFragger_3_Parser(
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
            "Raw data location": "/Users/cellzome/Dev/Gits/Ursgal/ursgal_master/example_data/test_Creinhardtii_QE_pH11.mzML",
            "15N": False,
        },
    )
    df = parser.unify()
    assert all(df["Modifications"] == "")


def test_map_mod_translation():
    assert 1 == 2


def testtranslate_mods():
    assert 1 == 2
