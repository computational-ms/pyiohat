#!/usr/bin/env python
import pandas as pd
import pytest
import xml.etree.ElementTree as etree

from pyprotista.parsers.ident.comet_2020_01_4_parser import (
    Comet_2020_01_4_Parser,
    get_modification_mass_mapping,
    get_xml_data,
    get_peptide_lookup,
    get_spec_records,
    map_mods_sequences,
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
    fixed_mods = {"C": "Carbamidomethyl"}
    modification_mass_map = {
        "57.021464": "Carbamidomethyl",
        "15.994915": "Oxidation",
        "42.010565": "Acetyl",
    }
    (version, spec_records,) = get_xml_data(
        file=input_file,
        mapping_dict=parser.mapping_dict,
        fixed_mods=fixed_mods,
        modification_mass_map=modification_mass_map,
    )

    assert version == "comet_2020_01_4"
    assert len(spec_records) == 60
    assert modification_mass_map["57.021464"] == "Carbamidomethyl"
    assert spec_records[24] == {
        "comet:num_matched_ions": "2",
        "comet:num_unmatched_ions": "22",
        "comet:xcorr": "0.2554",
        "comet:deltacn": "1.0000",
        "comet:score": "5.4000",
        "comet:e_value": "3.84E+01",
        "charge": "3",
        "sequence": "LCVLHEK",
        "exp_mz": "300.165802",
        "calc_mz": "300.165351",
        "modifications": "Carbamidomethyl:2",
        "spectrum_id": "2841",
    }
    assert spec_records[7] == {
        "comet:num_matched_ions": "1",
        "comet:num_unmatched_ions": "35",
        "comet:xcorr": "0.2159",
        "comet:deltacn": "1.0000",
        "comet:score": "2.8000",
        "comet:e_value": "2.42E+02",
        "charge": "3",
        "sequence": "ECCDKPLLEK",
        "exp_mz": "431.205261",
        "calc_mz": "431.205546",
        "modifications": "Carbamidomethyl:2;Carbamidomethyl:3",
        "spectrum_id": "2615",
    }
    df = pd.DataFrame(spec_records)
    unique_sequences = df["sequence"].unique().tolist()
    assert len(unique_sequences) == 21
    unique_mods = df["modifications"].unique().tolist()
    assert len(unique_mods) == 9
    assert "" in unique_mods
    assert "HLVDEPQNLIK" and "YICDNQDTISSK" and "LKPDPNTLCDEFK" in unique_sequences

    fixed_mods = {}
    modification_mass_map = {
        "57.021464": "Carbamidomethyl",
        "15.994915": "Oxidation",
        "42.010565": "Acetyl",
    }
    (version, spec_records,) = get_xml_data(
        file=input_file,
        mapping_dict=parser.mapping_dict,
        fixed_mods=fixed_mods,
        modification_mass_map=modification_mass_map,
    )
    assert version == "comet_2020_01_4"
    assert spec_records[1] == {
        "comet:num_matched_ions": "1",
        "comet:num_unmatched_ions": "13",
        "comet:xcorr": "0.0879",
        "comet:deltacn": "1.0000",
        "comet:score": "2.3000",
        "comet:e_value": "5.75E+02",
        "charge": "2",
        "sequence": "LRCASIQK",
        "exp_mz": "509.780029",
        "calc_mz": "509.279126",
        "modifications": "Acetyl:0",
        "spectrum_id": "2497",
    }
    df = pd.DataFrame(spec_records)
    unique_sequences = df["sequence"].unique().tolist()
    assert len(unique_sequences) == 21
    unique_mods = df["modifications"].unique().tolist()
    assert "Carbamidomethyl" not in unique_mods
    assert "Acetyl:0" in unique_mods
    assert len(unique_mods) == 2


def test_engine_parsers_comet_get_peptide_lookup():
    peptide = etree.Element("Peptide", id="LRCASIQK;8:42.010565;")
    peptide_sequence = etree.Element("PeptideSequence")
    peptide_sequence.text = "LRCASIQK"
    modification = etree.Element(
        "Modification", location="0", monoisotopicMassDelta="42.010565"
    )
    fixed_mods = {"C": "Carbamidomethyl"}
    modification_mass_map = {
        "57.021464": "Carbamidomethyl",
        "15.994915": "Oxidation",
        "42.010565": "Acetyl",
    }
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
            fixed_mods=fixed_mods,
            modification_mass_map=modification_mass_map,
        )

    assert len(peptide_lookup) == 1
    assert peptide_lookup["LRCASIQK;8:42.010565;"] == (
        "Carbamidomethyl:3;Acetyl:0",
        "LRCASIQK",
    )
    assert sequence == "LRCASIQK"
    assert modifications == "Carbamidomethyl:3;Acetyl:0"
    assert peptide_lookup == {
        "LRCASIQK;8:42.010565;": ("Carbamidomethyl:3;Acetyl:0", "LRCASIQK")
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
    peptide_lookup = {
        "AEFVEVTK;": ("", "AEFVEVTK"),
        "AWSVAR;": ("", "AWSVAR"),
        "CCTESLVNR;": ("Carbamidomethyl:1;Carbamidomethyl:2", "CCTESLVNR"),
        "DDSPDLPK;": ("", "DDSPDLPK"),
        "DLGEEHFK;": ("", "DLGEEHFK"),
        "DLGEEHFK;8:42.010565;": ("Acetyl:0", "DLGEEHFK"),
        "EACFAVEGPK;": ("Carbamidomethyl:3", "EACFAVEGPK"),
        "EACFAVEGPK;10:42.010565;": ("Carbamidomethyl:3;Acetyl:0", "EACFAVEGPK"),
        "ECCDKPLLEK;": ("Carbamidomethyl:2;Carbamidomethyl:3", "ECCDKPLLEK"),
        "GACLLPK;": ("Carbamidomethyl:3", "GACLLPK"),
        "HLVDEPQNLIK;": ("", "HLVDEPQNLIK"),
        "LCVLHEK;": ("Carbamidomethyl:2", "LCVLHEK"),
        "LGEYGFQNALIVR;": ("", "LGEYGFQNALIVR"),
        "LKPDPNTLCDEFK;": ("Carbamidomethyl:9", "LKPDPNTLCDEFK"),
        "LRCASIQK;8:42.010565;": ("Carbamidomethyl:3;Acetyl:0", "LRCASIQK"),
        "LVTDLTK;": ("", "LVTDLTK"),
        "LVTDLTK;7:42.010565;": ("Acetyl:0", "LVTDLTK"),
        "LVVSTQTALA;": ("", "LVVSTQTALA"),
        "QEPERNECFLSHK;": ("Carbamidomethyl:8", "QEPERNECFLSHK"),
        "SHCIAEVEK;": ("Carbamidomethyl:3", "SHCIAEVEK"),
        "VPQVSTPTLVEVSR;": ("", "VPQVSTPTLVEVSR"),
        "YICDNQDTISSK;": ("Carbamidomethyl:3", "YICDNQDTISSK"),
        "YLYEIAR;": ("", "YLYEIAR"),
        "YLYEIARR;": ("", "YLYEIARR"),
    }

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
            peptide_lookup=peptide_lookup,
        )
    assert spec_records == [
        {
            "comet:num_matched_ions": "3",
            "charge": "3",
            "sequence": "SHCIAEVEK",
            "exp_mz": "358.174682",
            "calc_mz": "358.174575",
            "modifications": "Carbamidomethyl:3",
            "spectrum_id": "2458",
        }
    ]
    assert len(spec_records[0]) == 7
    assert len(spec_records) == 1


def test_engine_parsers_comet_get_modification_mass_mapping():
    input_file = pytest._test_path / "data" / "BSA1_comet_2020_01_4.mzid"
    fixed_mods, modification_mass_map = get_modification_mass_mapping(input_file)
    assert fixed_mods == {"C": "Carbamidomethyl"}
    assert len(fixed_mods) == 1
    assert modification_mass_map == {
        "57.021464": "Carbamidomethyl",
        "15.994915": "Oxidation",
        "42.010565": "Acetyl",
    }
    assert len(modification_mass_map) == 3


def test_engine_parsers_comet_map_mods_sequences():
    fixed_mods = {"C": "Carbamidomethyl", "O": "Oxidation"}
    sequence = "ABCDEFCABC"
    mods = map_mods_sequences(fixed_mods, sequence)
    assert mods == "Carbamidomethyl:3;Carbamidomethyl:7;Carbamidomethyl:10"
    sequence = "ABCDOEFOCABC"
    mods = map_mods_sequences(fixed_mods, sequence)
    assert (
        mods
        == "Carbamidomethyl:3;Carbamidomethyl:9;Carbamidomethyl:12;Oxidation:5;Oxidation:8"
    )
    sequence = "OOOOOCCCCCMMMMM"
    mods = map_mods_sequences(fixed_mods, sequence)
    assert (
        mods
        == "Carbamidomethyl:6;Carbamidomethyl:7;Carbamidomethyl:8;Carbamidomethyl:9;Carbamidomethyl:10;Oxidation:1;Oxidation:2;Oxidation:3;Oxidation:4;Oxidation:5"
    )
    sequence = "CD"
    mods = map_mods_sequences(fixed_mods, sequence)
    assert mods == "Carbamidomethyl:1"
    sequence = "C"
    mods = map_mods_sequences(fixed_mods, sequence)
    assert mods == "Carbamidomethyl:1"
    fixed_mods = {}
    sequence = "ASLDPOCSADK"
    with pytest.raises(ValueError):
        map_mods_sequences(fixed_mods, sequence)
