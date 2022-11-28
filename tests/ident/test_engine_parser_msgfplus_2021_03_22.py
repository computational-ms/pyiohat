#!/usr/bin/env python

import pytest

from pyprotista.parsers.ident.msgfplus_2021_03_22_parser import (
    MSGFPlus_2021_03_22_Parser,
    get_version,
    get_peptide_lookup,
    get_spec_records,
)


def test_engine_parsers_msgfplus_init():
    input_file = pytest._test_path / "data" / "BSA1_msgfplus_2021_03_22.mzid"
    parser = MSGFPlus_2021_03_22_Parser(
        input_file,
        params={
            "cpus": 2,
            "enzyme": "(?<=[KR])(?![P])",
            "terminal_cleavage_site_integrity": "any",
            "validation_score_field": {"msgfplus_2021_03_22": "ms-gf:spec_evalue"},
            "bigger_scores_better": {"msgfplus_2021_03_22": False},
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


def test_engine_parsers_msgfplus_check_parser_compatibility():
    msgf_parser_class = MSGFPlus_2021_03_22_Parser
    input_file = pytest._test_path / "data" / "BSA1_msgfplus_2021_03_22.mzid"
    assert msgf_parser_class.check_parser_compatibility(input_file) is True


def test_engine_parsers_msgfplus_check_parser_compatibility_fail_with_omssa_file():
    msgf_parser_class = MSGFPlus_2021_03_22_Parser
    input_file = (
        pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_omssa_2_1_9.csv"
    )
    assert msgf_parser_class.check_parser_compatibility(input_file) is False


def test_engine_parsers_msgfplus_check_dataframe_integrity():
    input_file = pytest._test_path / "data" / "BSA1_msgfplus_2021_03_22.mzid"
    rt_lookup_path = pytest._test_path / "data" / "BSA1_ursgal_lookup.csv"
    db_path = pytest._test_path / "data" / "BSA.fasta"

    parser = MSGFPlus_2021_03_22_Parser(
        input_file,
        params={
            "cpus": 2,
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "enzyme": "(?<=[KR])(?![P])",
            "terminal_cleavage_site_integrity": "any",
            "validation_score_field": {"msgfplus_2021_03_22": "ms-gf:spec_evalue"},
            "bigger_scores_better": {"msgfplus_2021_03_22": False},
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
    print(df)
    assert pytest.approx(df["exp_mz"].mean()) == 488.0319
    assert len(df) == 92
    assert pytest.approx(df["ucalc_mz"].mean()) == 488.03167
    assert (df["raw_data_location"] == "path/for/glory.mzML").all()
    assert df["modifications"].str.contains("Acetyl:0").sum() == 0
    assert df["modifications"].str.contains("Oxidation:").sum() == 0
    assert (
        df["modifications"].str.count("Carbamidomethyl:")
        == df["sequence"].str.count("C")
    ).all()
    assert df["modifications"].str.count(":").sum() == 71


def test_engine_parsers_msgfplus_check_dataframe_integrity_unknown_mod():
    input_file = (
        pytest._test_path / "data" / "BSA1_msgfplus_2021_03_22_unknown_mod.mzid"
    )
    rt_lookup_path = pytest._test_path / "data" / "BSA1_ursgal_lookup.csv"
    db_path = pytest._test_path / "data" / "BSA.fasta"

    parser = MSGFPlus_2021_03_22_Parser(
        input_file,
        params={
            "cpus": 2,
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "enzyme": "(?<=[KR])(?![P])",
            "terminal_cleavage_site_integrity": "any",
            "validation_score_field": {"msgfplus_2021_03_22": "ms-gf:spec_evalue"},
            "bigger_scores_better": {"msgfplus_2021_03_22": False},
            "modifications": [
                {
                    "aa": "M",
                    "type": "opt",
                    "position": "any",
                    "name": "Oxidation",
                },
                {
                    "aa": "C",
                    "type": "opt",
                    "position": "any",
                    "name": "Carbamidomethyl",
                },
                {
                    "aa": "C",
                    "type": "opt",
                    "position": "any",
                    "name": "DTB-IAA",
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

    assert len(df) == 92
    # Currently this excludes two peptides and is reduced
    assert pytest.approx(df["ucalc_mz"].mean()) == 486.56
    assert (df["raw_data_location"] == "path/for/glory.mzML").all()
    assert df["modifications"].str.contains("Acetyl:0").sum() == 0
    assert df["modifications"].str.contains("Oxidation:").sum() == 0
    # the unknown modification is actually DTB-IAA
    assert df["modifications"].str.contains("DTB-IAA:3").sum() == 2
    assert (
        df[df["sequence"] != "EACFAVEGPK"]["modifications"].str.count(
            "Carbamidomethyl:"
        )
        == df[df["sequence"] != "EACFAVEGPK"]["sequence"].str.count("C")
    ).all()
    assert df["modifications"].str.count(":").sum() == 71


def test_engine_parsers_msgfplus_get_version():
    input_file = pytest._test_path / "data" / "BSA1_msgfplus_2021_03_22.mzid"
    version = get_version(input_file)

    assert version == "msgfplus_2021_03_22"


def test_engine_parsers_msgfplus_get_peptide_lookup():
    input_file = pytest._test_path / "data" / "BSA1_msgfplus_2021_03_22.mzid"
    peptide_lookup = get_peptide_lookup(input_file)

    assert len(peptide_lookup) == 24
    for pep in peptide_lookup:
        assert len(peptide_lookup[pep]) == 2
        assert pep.startswith("Pep_")
    assert peptide_lookup["Pep_LVTDLTK"] == {"sequence": "LVTDLTK", "modifications": ""}
    assert peptide_lookup["Pep_LVVSTQTALA"] == {
        "sequence": "LVVSTQTALA",
        "modifications": "",
    }


def test_engine_parsers_msgfplus_get_spec_records():
    input_file = pytest._test_path / "data" / "BSA1_msgfplus_2021_03_22.mzid"
    obj = MSGFPlus_2021_03_22_Parser(input_file=None, params=None)
    peptide_lookup = {
        "Pep_YICDNQDTISSK": {
            "sequence": "YICDNQDTISSK",
            "modifications": "Carbamidomethyl:3",
        },
        "Pep_CCTESLVNR": {
            "sequence": "CCTESLVNR",
            "modifications": "Carbamidomethyl:1;Carbamidomethyl:2",
        },
        "Pep_EYEATLEECCAK": {
            "sequence": "EYEATLEECCAK",
            "modifications": "Carbamidomethyl:9;Carbamidomethyl:10",
        },
        "Pep_ETYGDMADCCEK": {
            "sequence": "ETYGDMADCCEK",
            "modifications": "Carbamidomethyl:9;Carbamidomethyl:10",
        },
        "Pep_EACFAVEGPK": {
            "sequence": "EACFAVEGPK",
            "modifications": "Carbamidomethyl:3",
        },
        "Pep_DDSPDLPK": {"sequence": "DDSPDLPK", "modifications": ""},
        "Pep_ECCDKPLLEK": {
            "sequence": "ECCDKPLLEK",
            "modifications": "Carbamidomethyl:2;Carbamidomethyl:3",
        },
        "Pep_HLVDEPQNLIK": {"sequence": "HLVDEPQNLIK", "modifications": ""},
        "Pep_DLGEEHFK": {"sequence": "DLGEEHFK", "modifications": ""},
        "Pep_DDPHACYSTVFDK": {
            "sequence": "DDPHACYSTVFDK",
            "modifications": "Carbamidomethyl:6",
        },
        "Pep_AEFVEVTK": {"sequence": "AEFVEVTK", "modifications": ""},
        "Pep_YLYEIAR": {"sequence": "YLYEIAR", "modifications": ""},
        "Pep_LCVLHEK": {"sequence": "LCVLHEK", "modifications": "Carbamidomethyl:2"},
        "Pep_LVTDLTK": {"sequence": "LVTDLTK", "modifications": ""},
        "Pep_YNGVFQECCQAEDK": {
            "sequence": "YNGVFQECCQAEDK",
            "modifications": "Carbamidomethyl:8;Carbamidomethyl:9",
        },
        "Pep_LKPDPNTLCDEFK": {
            "sequence": "LKPDPNTLCDEFK",
            "modifications": "Carbamidomethyl:9",
        },
        "Pep_SHCIAEVEK": {
            "sequence": "SHCIAEVEK",
            "modifications": "Carbamidomethyl:3",
        },
        "Pep_QEPERNECFLSHK": {
            "sequence": "QEPERNECFLSHK",
            "modifications": "Carbamidomethyl:8",
        },
        "Pep_VPQVSTPTLVEVSR": {"sequence": "VPQVSTPTLVEVSR", "modifications": ""},
        "Pep_LVVSTQTALA": {"sequence": "LVVSTQTALA", "modifications": ""},
        "Pep_AWSVAR": {"sequence": "AWSVAR", "modifications": ""},
        "Pep_GACLLPK": {"sequence": "GACLLPK", "modifications": "Carbamidomethyl:3"},
        "Pep_YLYEIARR": {"sequence": "YLYEIARR", "modifications": ""},
        "Pep_LGEYGFQNALIVR": {"sequence": "LGEYGFQNALIVR", "modifications": ""},
    }

    spec_records = get_spec_records(input_file, peptide_lookup, obj.mapping_dict)

    assert len(spec_records) == 92
    assert spec_records[0] == {
        "ms-gf:raw_score": "40",
        "ms-gf:denovoscore": "40",
        "ms-gf:spec_evalue": "4.4458354E-15",
        "ms-gf:evalue": "2.6986221E-12",
        "ms-gf:num_matched_ions": "3",
        "charge": "2",
        "exp_mz": "722.3272094726562",
        "calc_mz": "722.3246459960938",
        "sequence": "YICDNQDTISSK",
        "modifications": "Carbamidomethyl:3",
        "spectrum_title": "glory.2791.2791.2",
        "spectrum_id": "2791",
        "retention_time_seconds": "1918.6086",
    }
    assert spec_records[44] == {
        "ms-gf:raw_score": "26",
        "ms-gf:denovoscore": "32",
        "ms-gf:spec_evalue": "2.2523169E-8",
        "ms-gf:evalue": "1.36715635E-5",
        "ms-gf:num_matched_ions": "1",
        "charge": "3",
        "exp_mz": "325.4911804199219",
        "calc_mz": "325.4907531738281",
        "sequence": "DLGEEHFK",
        "modifications": "",
        "spectrum_title": "glory.mzML.2663.2663.3",
        "spectrum_id": "2663",
        "retention_time_seconds": "1840.7933333333335",
    }
