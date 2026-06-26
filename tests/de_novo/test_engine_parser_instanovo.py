import pytest

from pyiohat.parsers.de_novo.instanovo_1_parser import Instanovo_1_Parser


def test_engine_parsers_instanovo_init():
    input_file = pytest._test_path / "data" / "test_instanovo.instanovo.csv"
    parser = Instanovo_1_Parser(
        input_file,
        params={
            "cpus": 2,
            "enzyme": "(?<=[KR])(?![P])",
            "terminal_cleavage_site_integrity": "any",
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


def test_engine_parsers_instanovo_metadata():
    input_file = pytest._test_path / "data" / "test_instanovo.instanovo.csv"
    parser = Instanovo_1_Parser(
        input_file,
        params={
            "cpus": 2,
            "enzyme": "(?<=[KR])(?![P])",
            "terminal_cleavage_site_integrity": "any",
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
    assert parser.metadata


def test_engine_parsers_instanovo_check_parser_compatibility():
    input_file = pytest._test_path / "data" / "test_instanovo.instanovo.csv"
    assert Instanovo_1_Parser.check_parser_compatibility(input_file) is True


def test_engine_parsers_instanovo_check_dataframe_integrity():
    input_file = pytest._test_path / "data" / "test_instanovo.instanovo.csv"
    rt_lookup_path = pytest._test_path / "data" / "instanovo_lookup.csv"
    db_path = pytest._test_path / "data" / "Hfvol_prot_250410.fasta"

    parser = Instanovo_1_Parser(
        input_file,
        params={
            "cpus": 2,
            "database": db_path,
            "rt_pickle_name": rt_lookup_path,
            "enzyme": "(?<=[KR])(?!P)",
            "terminal_cleavage_site_integrity": "any",
            "min_pep_length": 5,
            "modifications": [
                {
                    "aa": "M",
                    "type": "opt",
                    "position": "any",
                    "name": "Oxidation",
                },
                {
                    "aa": "S",
                    "type": "opt",
                    "position": "any",
                    "name": "Phospho",
                },
                {
                    "aa": "T",
                    "type": "opt",
                    "position": "any",
                    "name": "Phospho",
                },
                {
                    "aa": "Y",
                    "type": "opt",
                    "position": "any",
                    "name": "Phospho",
                },
                {
                    "aa": "Q",
                    "type": "opt",
                    "position": "any",
                    "name": "Deamidated",
                },
                {
                    "aa": "N",
                    "type": "opt",
                    "position": "any",
                    "name": "Deamidated",
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
                {
                    "aa": "*",
                    "type": "opt",
                    "position": "any",
                    "name": "Carbamyl",
                },
                {
                    "aa": "*",
                    "type": "opt",
                    "position": "any",
                    "name": "Ammonia-loss",
                },
            ],
        },
    )
    df = parser.unify()
    assert len(df) == 100
    assert pytest.approx(df["ucalc_mz"].mean(), abs=1e-3) == 370.133
    assert pytest.approx(df["exp_mz"].mean(), abs=1e-3) == 471.923

    assert df["modifications"].str.contains("Carbamidomethyl:2").sum() == 1
    assert df["modifications"].str.contains("Oxidation:").sum() == 17
    assert (
        df["modifications"].str.count("Carbamidomethyl:")
        == df["sequence"].str.count("C")
    ).all()
    assert df["modifications"].str.count(":").sum() == 23
    assert (df["raw_data_location"] == "path/for/glory.mzML").all()


def test_map_mod_translation_instanovo():
    import pandas as pd

    input_file = pytest._test_path / "data" / "test_instanovo.instanovo.csv"
    db_path = pytest._test_path / "data" / "Hfvol_prot_250410.fasta"
    rt_lookup_path = pytest._test_path / "data" / "instanovo_lookup.csv"

    parser = Instanovo_1_Parser(
        input_file,
        params={
            "cpus": 2,
            "database": db_path,
            "rt_pickle_name": rt_lookup_path,
            "enzyme": "(?<=[KR])(?!P)",
            "terminal_cleavage_site_integrity": "any",
            "min_pep_length": 5,
            "modifications": [
                {
                    "aa": "M",
                    "type": "opt",
                    "position": "any",
                    "name": "Oxidation",
                },
                {
                    "aa": "S",
                    "type": "opt",
                    "position": "any",
                    "name": "Phospho",
                },
                {
                    "aa": "T",
                    "type": "opt",
                    "position": "any",
                    "name": "Phospho",
                },
                {
                    "aa": "Y",
                    "type": "opt",
                    "position": "any",
                    "name": "Phospho",
                },
                {
                    "aa": "Q",
                    "type": "opt",
                    "position": "any",
                    "name": "Deamidated",
                },
                {
                    "aa": "N",
                    "type": "opt",
                    "position": "any",
                    "name": "Deamidated",
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
                {
                    "aa": "*",
                    "type": "opt",
                    "position": "any",
                    "name": "Carbamyl",
                },
                {
                    "aa": "*",
                    "type": "opt",
                    "position": "any",
                    "name": "Ammonia-loss",
                },
            ],
        },
    )
    # Manual injection of a dataframe to test specific mapping logic
    parser.df.loc[5, "sequence"] = "M[UNIMOD:35]FEKKFK"
    df = parser.unify()

    assert df["modifications"].iloc[5] == "Oxidation:1"
    assert df["sequence"].iloc[5] == "MFEKKFK"


def test_instanovo_parser_missing_ms_level_error():
    input_file = pytest._test_path / "data" / "test_instanovo.instanovo.csv"
    db_path = pytest._test_path / "data" / "Hfvol_prot_250410.fasta"
    bad_lookup_path = pytest._test_path / "data" / "instanovo_lookup_no_ms_level.csv"
    with pytest.raises(KeyError) as excinfo:
        parser = Instanovo_1_Parser(
            input_file,
            params={
                "cpus": 2,
                "database": db_path,
                "rt_pickle_name": bad_lookup_path,
                "enzyme": "(?<=[KR])(?!P)",
                "modifications": [],
                "terminal_cleavage_site_integrity": "any",
                "min_pep_length": 5,
            },
        )
        parser.unify()
    assert "Could not uniquely assign meta data" in str(excinfo.value)


def test_map_multiple_mod_translations_instanovo():
    """Test mapping and 1-based positioning when a single peptide sequence

    contains multiple, differing modifications (e.g. N-terminal and internal).
    """
    import pandas as pd

    input_file = pytest._test_path / "data" / "test_instanovo.instanovo.csv"
    db_path = pytest._test_path / "data" / "Hfvol_prot_250410.fasta"
    rt_lookup_path = pytest._test_path / "data" / "instanovo_lookup.csv"

    parser = Instanovo_1_Parser(
        input_file,
        params={
            "cpus": 2,
            "database": db_path,
            "rt_pickle_name": rt_lookup_path,
            "enzyme": "(?<=[KR])(?!P)",
            "terminal_cleavage_site_integrity": "any",
            "min_pep_length": 5,
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
    # Inject a mock row at index 0 containing multiple modifications:
    parser.df.loc[0, "sequence"] = "[UNIMOD:1]-M[UNIMOD:35]ACG[UNIMOD:4]K"

    df = parser.unify()
    # Get the processed row
    processed_row = df.iloc[0]
    # Check that the raw bracket codes are cleanly stripped out of the sequence
    assert processed_row["sequence"] == "MACGK"
    # Validate that multiple modifications are mapped to their correct 1-based string indices
    mod_string = processed_row["modifications"]
    assert "Acetyl:0" in mod_string
    assert "Oxidation:1" in mod_string
    assert "Carbamidomethyl:4" in mod_string


def test_instanovo_parser_mod_mapper_empty_list_raises_keyerror():
    """If mod_mapper.id_to_name() returns an empty list (lookup succeeded but
    found no matching name), parse_sequence_modifications should raise KeyError
    rather than silently dropping the modification.
    """
    from unittest.mock import patch
    import pandas as pd

    input_file = pytest._test_path / "data" / "test_instanovo.instanovo.csv"
    db_path = pytest._test_path / "data" / "Hfvol_prot_250410.fasta"
    rt_lookup_path = pytest._test_path / "data" / "instanovo_lookup.csv"

    parser = Instanovo_1_Parser(
        input_file,
        params={
            "cpus": 2,
            "database": db_path,
            "rt_pickle_name": rt_lookup_path,
            "enzyme": "(?<=[KR])(?!P)",
            "terminal_cleavage_site_integrity": "any",
            "min_pep_length": 5,
            "modifications": [
                {
                    "aa": "M",
                    "type": "opt",
                    "position": "any",
                    "name": "Oxidation",
                },
            ],
        },
    )

    parser.df.loc[0, "sequence"] = "M[UNIMOD:35]FEKKFK"

    with patch.object(parser.mod_mapper, "id_to_name", return_value=[]):
        with pytest.raises(KeyError, match="Unable to map"):
            parser.parse_sequence_modifications()
