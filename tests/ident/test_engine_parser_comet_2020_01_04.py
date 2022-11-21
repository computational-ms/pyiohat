#!/usr/bin/env python

import pytest
import xml.etree.ElementTree as etree

from pyprotista.parsers.ident.comet_2020_01_4_parser import (
    Comet_2020_01_4_Parser,
    get_xml_data,
    get_spec_records,
    get_peptide_lookup,
    get_modification_mass_map,
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


def test_engine_parsers_comet_get_xml_data():
    input_file = pytest._test_path / "data" / "BSA1_comet_2020_01_4.mzid"
    parser = Comet_2020_01_4_Parser(input_file=None, params=None)
    (
        version,
        peptide_lookup,
        spec_records,
        modification_mass_map,
        fixed_mods,
    ) = get_xml_data(input_file, parser.mapping_dict)

    assert len(peptide_lookup) == 24
    assert len(modification_mass_map) == 3
    assert len(spec_records) == 60
    assert modification_mass_map["57.021464"] == "Carbamidomethyl"
    assert (
        peptide_lookup["EACFAVEGPK;10:42.010565;"]["modifications"][0][
            "monoisotopicMassDelta"
        ]
        == "42.010565"
    )
    assert spec_records[0]["comet:e_value"] == "3.76E+01"


def test_engine_parsers_comet_get_peptide_lookup():
    peptide = etree.Element("Peptide", id="LRCASIQK;8:42.010565;")
    peptide_sequence = etree.Element("PeptideSequence")
    peptide_sequence.text = "LRCASIQK"
    modification = etree.Element(
        "Modification", location="0", monoisotopicMassDelta="42.010565"
    )
    results = [peptide_sequence, modification, peptide]

    modifications = []
    sequence = {}
    peptide_lookup = {}

    for entry in results:
        entry_tag = entry.tag
        sequence, modifications, peptide_lookup = get_peptide_lookup(
            entry=entry,
            entry_tag=entry_tag,
            sequence=sequence,
            modifications=modifications,
            peptide_lookup=peptide_lookup,
        )

    assert len(peptide_lookup) == 1
    assert peptide_lookup["LRCASIQK;8:42.010565;"]["modifications"][0] == {
        "monoisotopicMassDelta": "42.010565",
        "location": "0",
    }
    assert peptide_lookup["LRCASIQK;8:42.010565;"]["sequence"] == "LRCASIQK"
    assert peptide_lookup == {
        "LRCASIQK;8:42.010565;": {
            "modifications": [{"monoisotopicMassDelta": "42.010565", "location": "0"}],
            "sequence": "LRCASIQK",
        }
    }


def test_engine_parsers_comet_get_spec_records():
    cv_param = etree.Element(
        "cvParam",
        cvRef="MS",
        accession="MS:1001121",
        name="number of matched peaks",
        value="3",
    )
    spectrum_identification_item = etree.Element(
        "SpectrumIdentificationItem",
        id="SII_0",
        rank="0",
        chargeState="3",
        peptide_ref="SHCIAEVEK;",
        experimentalMassToCharge="358.174682",
        calculatedMassToCharge="358.174575",
        passThreshold="true",
    )
    spectrum_identification_result = etree.Element(
        "SpectrumIdentificationResult",
        id="SIR_16",
        spectrumID="scan=2458",
        spectraData_ref="SD0",
    )
    results = [cv_param, spectrum_identification_item, spectrum_identification_result]
    parser = Comet_2020_01_4_Parser(input_file=None, params=None)

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
            mapping_dict=parser.mapping_dict,
        )
    assert spec_records == [
        {
            "comet:num_matched_ions": "3",
            "charge": "3",
            "sequence": "SHCIAEVEK;",
            "exp_mz": "358.174682",
            "calc_mz": "358.174575",
            "spectrum_id": "2458",
        }
    ]
    assert len(spec_records[0]) == 6
    assert spec_records[0]["comet:num_matched_ions"] == "3"
    assert spec_records[0]["charge"] == "3"
    assert spec_records[0]["sequence"] == "SHCIAEVEK;"
    assert spec_records[0]["exp_mz"] == "358.174682"
    assert spec_records[0]["calc_mz"] == "358.174575"
    assert spec_records[0]["spectrum_id"] == "2458"


def test_engine_parsers_comet_get_modification_mass_map():
    cv_param = etree.Element(
        "cvParam", cvRef="UNIMOD", accession="UNIMOD:UNIMOD:4", name="Carbamidomethyl"
    )
    search_modification = etree.Element(
        "SearchModification", residues="C", massDelta="57.021464", fixedMod="true"
    )
    results = [cv_param, search_modification]

    modification_mass_map = {}
    mod_name = ""
    fixed_mods = {}

    for entry in results:
        entry_tag = entry.tag
        modification_mass_map, mod_name, fixed_mods = get_modification_mass_map(
            entry=entry,
            entry_tag=entry_tag,
            modification_mass_map=modification_mass_map,
            mod_name=mod_name,
            fixed_mods=fixed_mods,
        )

    assert fixed_mods == {"C": "Carbamidomethyl"}
    assert mod_name == "Carbamidomethyl"
    assert modification_mass_map == {"57.021464": "Carbamidomethyl"}
