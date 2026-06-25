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
    # 1. Acetyl at the N-terminus ([UNIMOD:1]-...)
    # 2. Oxidation at the first Methionine (index 1)
    # 3. Carbamidomethylation at Cysteine (index 4)
    # The clean sequence should be: MACGKR
    parser.df.loc[0, "sequence"] = "[UNIMOD:1]-M[UNIMOD:35]ACG[UNIMOD:4]K"

    df = parser.unify()

    # Get our processed row
    processed_row = df.iloc[0]

    # 1. Check that the raw bracket codes are cleanly stripped out of the sequence
    assert processed_row["sequence"] == "MACGK"

    # 2. Validate that multiple modifications are mapped to their correct 1-based string indices
    # Expecting: Acetyl at position 0 (or 1 depending on N-term style), Oxidation at M (1), Carbamidomethyl at C (3)
    # Adjust string formatting expectations based on how your base parser joins multiple mods (usually semi-colon separated)
    mod_string = processed_row["modifications"]
    assert "Acetyl:0" in mod_string
    assert "Oxidation:1" in mod_string
    assert "Carbamidomethyl:4" in mod_string


# Tests for n-term and digits
# def test_c_terminal_tmt():
#     input_file = pytest._test_path / "data" / "test_instanovo.instanovo.csv"

#     parser = Instanovo_1_Parser(
#         input_file,
#         params={
#             "cpus": 2,
#             "enzyme": "(?<=[KR])(?![P])",
#             "terminal_cleavage_site_integrity": "any",
#             "modifications": [
#                 {
#                     "aa": "M",
#                     "type": "opt",
#                     "position": "any",
#                     "name": "Oxidation",
#                 },
#                 {
#                     "aa": "*",
#                     "type": "opt",
#                     "position": "Prot-N-term",
#                     "name": "Acetyl",
#                 },
#                 {
#                     "aa": "*",
#                     "type": "opt",
#                     "position": "N-term",
#                     "name": "TMT6plex",
#                 },
#                 {
#                     "aa": "C",
#                     "type": "fix",
#                     "position": "any",
#                     "name": "Carbamidomethyl",
#                 },
#                 {
#                     "aa": "K",
#                     "type": "fix",
#                     "position": "any",
#                     "name": "TMT6plex",
#                 },
#             ],
#             "15N": False,
#         },
#     )
#     converted = parser.translate_mods()
#     assert converted[0] == "TMT6plex:0;TMT6plex:6"
#     assert converted[1] == "TMT6plex:0"


# def test_msfragger_open_search():
#    input_file = pytest._test_path / "data" / "BSA1_open_search.msfragger4.tsv"
#    rt_lookup_path = pytest._test_path / "data" / "BSA1_ursgal_lookup.csv"
#    db_path = pytest._test_path / "data" / "BSA.fasta"

#    parser = Instanovo_1_Parser(
#        input_file,
#        params={
#            "cpus": 2,
#            "rt_pickle_name": rt_lookup_path,
#            "database": db_path,
#            "enzyme": "(?<=[KR])(?![P])",
#            "terminal_cleavage_site_integrity": "any",
#            "validation_score_field": {"MSFragger_4_0": "msfragger:hyperscore"},
#            "bigger_scores_better": {"MSFragger_4_0": True},
#            "modifications": [],
#            "15N": False,
#        },
#    )
#    df = parser.unify()
#    assert df["mass_delta"].mean() == pytest.approx(458.901, abs=1e-6)
