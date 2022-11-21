#!/usr/bin/env python

import pytest
import xml.etree.ElementTree as etree

from pyprotista.parsers.ident.msgfplus_2021_03_22_parser import (
    MSGFPlus_2021_03_22_Parser,
    get_xml_data,
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


def test_engine_parsers_msgfplus_get_xml_data():
    input_file = pytest._test_path / "data" / "BSA1_msgfplus_2021_03_22.mzid"
    obj = MSGFPlus_2021_03_22_Parser(input_file=None, params=None)

    version, peptide_lookup, spec_records = get_xml_data(input_file, obj.mapping_dict)

    assert version == "msgfplus_2021_03_22"
    assert len(spec_records) == 92
    assert len(peptide_lookup) == 24
    assert spec_records[3]["ms-gf:evalue"] == "1.6910947E-9"
    assert spec_records[3]["retention_time_seconds"] == "1793.826"
    assert (
        peptide_lookup["Pep_CCTESLVNR"]["modifications"]
        == "Carbamidomethyl:1;Carbamidomethyl:2"
    )


def test_engine_parsers_msgfplus_get_peptide_lookup():

    cv_param = etree.Element(
        "cvParam", cvRef="UNIMOD", accession="UNIMOD:4", name="Carbamidomethyl"
    )
    modification = etree.Element(
        "Modification", location="6", monoisotopicMassDelta="57.021464"
    )
    peptide_sequence = etree.Element("PeptideSequence")
    peptide_sequence.text = "DDPHACYSTVFDK"
    peptide = etree.Element("Peptide", id="Pep_DDPHACYSTVFDK")

    results = [cv_param, modification, peptide_sequence, peptide]

    cv_param_modifications = ""
    sequence = {}
    peptide_lookup = {}

    for entry in results:
        entry_tag = entry.tag
        sequence, cv_param_modifications, peptide_lookup = get_peptide_lookup(
            entry=entry,
            entry_tag=entry_tag,
            sequence=sequence,
            cv_param_modifications=cv_param_modifications,
            peptide_lookup=peptide_lookup,
        )

    assert len(peptide_lookup) == 1
    assert peptide_lookup["Pep_DDPHACYSTVFDK"]["modifications"] == "Carbamidomethyl:6"
    assert peptide_lookup["Pep_DDPHACYSTVFDK"]["sequence"] == "DDPHACYSTVFDK"
    assert cv_param_modifications == ""


def test_engine_parsers_msgfplus_get_spec_records():
    cv_param = etree.Element(
        "cvParam",
        cvRef="PSI-MS",
        accession="MS:1002049",
        name="MS-GF:RawScore",
        value="40",
    )
    user_param = etree.Element("userParam", name="NumMatchedMainIons", value="3")
    spectrum_identification_item = etree.Element(
        "SpectrumIdentificationItem",
        chargeState="2",
        experimentalMassToCharge="722.3272094726562",
        calculatedMassToCharge="722.3246459960938",
        peptide_ref="Pep_YICDNQDTISSK",
        rank="1",
        passThreshold="true",
        id="SII_350_1",
    )
    spectrum_identification_result = etree.Element(
        "SpectrumIdentificationResult",
        spectrumID="index=349",
        spectraData_ref="SID_1",
        id="SIR_350",
    )

    results = [
        cv_param,
        user_param,
        spectrum_identification_item,
        spectrum_identification_result,
    ]

    obj = MSGFPlus_2021_03_22_Parser(input_file=None, params=None)
    spec_results = {}
    spec_ident_items = []
    spec_records = []

    for entry in results:
        entry_tag = entry.tag
        spec_results, spec_ident_items, spec_records = get_spec_records(
            entry=entry,
            entry_tag=entry_tag,
            spec_results=spec_results,
            spec_ident_items=spec_ident_items,
            spec_records=spec_records,
            mapping_dict=obj.mapping_dict,
        )

    assert len(spec_records) == 1
    assert len(spec_records[0]) == 6
    assert len(spec_ident_items) == 0
    assert spec_records[0]["exp_mz"] == "722.3272094726562"
