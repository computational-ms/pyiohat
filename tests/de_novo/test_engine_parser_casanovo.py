import pytest

from pyiohat.parsers.de_novo.casanovo_5_parser import Casanovo_5_Parser


def test_engine_parsers_casanovo_init():
    input_file = pytest._test_path / "data" / "test_casanovo.mztab"
    parser = Casanovo_5_Parser(
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


def test_engine_parsers_casanovo_metadata():
    input_file = pytest._test_path / "data" / "test_casanovo.mztab"
    parser = Casanovo_5_Parser(
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


def test_engine_parsers_casanovo_check_parser_compatibility():
    input_file = pytest._test_path / "data" / "test_casanovo.mztab"
    assert Casanovo_5_Parser.check_parser_compatibility(input_file) is True


def test_engine_parsers_casanovo_check_dataframe_integrity():
    input_file = pytest._test_path / "data" / "test_casanovo.mztab"
    rt_lookup_path = pytest._test_path / "data" / "casanovo_lookup.csv"
    db_path = pytest._test_path / "data" / "Hfvol_prot_250410.fasta"

    parser = Casanovo_5_Parser(
        input_file,
        params={
            "cpus": 2,
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
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
            "15N": False,
        },
    )
    df = parser.unify()
    assert len(df) == 101
    assert pytest.approx(df["ucalc_mz"].mean(), abs=1e-3) == 456.160
    assert pytest.approx(df["exp_mz"].mean(), abs=1e-3) == 456.334

    assert df["modifications"].str.contains("Carbamyl:0").sum() == 1
    assert df["modifications"].str.contains("Oxidation:").sum() == 12
    assert (
        df["modifications"].str.count("Carbamidomethyl:")
        == df["sequence"].str.count("C")
    ).all()
    assert df["modifications"].str.count(":").sum() == 21
    assert (df["raw_data_location"] == "path/for/glory.mzML").all()


def test_map_mod_translation_casanovo():
    import pandas as pd

    input_file = pytest._test_path / "data" / "test_casanovo.mztab"
    db_path = pytest._test_path / "data" / "Hfvol_prot_250410.fasta"
    rt_lookup_path = pytest._test_path / "data" / "casanovo_lookup.csv"

    parser = Casanovo_5_Parser(
        input_file,
        params={
            "cpus": 2,
            "database": db_path,
            "rt_pickle_name": rt_lookup_path,
            "enzyme": "(?<=[KR])(?!P)",
            "terminal_cleavage_site_integrity": "any",
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

    # Manual injection of a dataframe to test specific mapping logic
    parser.df = pd.DataFrame(
        {
            "sequence": ["LM[Oxidation]DKPEQLR"],
            "spectrum_id": ["ms_run[1]:controllerType=0 controllerNumber=1 scan=917"],
            "charge": [2],
            "casanovo:search_engine_score[1]": [0.99],
            "protein_id": [None],
            "exp_mz": [447.34625],
            "retention_time_seconds": [213.85],
        }
    )
    df = parser.unify()

    assert df["modifications"].iloc[0] == "Oxidation:2"
    assert df["sequence"].iloc[0] == "LMDKPEQLR"


# Tests for n-term and digits
# def test_c_terminal_tmt():
#     input_file = pytest._test_path / "data" / "test_casanovo.mztab"

#     parser = Casanovo_5_Parser(
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

#    parser = Casanovo_5_Parser(
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
