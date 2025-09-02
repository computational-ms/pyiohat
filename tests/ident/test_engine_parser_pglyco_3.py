import pytest

from pyiohat.parsers.ident.pglyco_3_parser import PGlyco_3_Parser


def test_engine_parsers_pglyco_init():
    input_file = pytest._test_path / "data" / "test_pglyco_3.txt"
    parser = PGlyco_3_Parser(
        input_file,
        params={
            "cpus": 2,
            "enzyme": "(?<=[KR])(?![P])",
            "terminal_cleavage_site_integrity": "any",
            # "validation_score_field": {"MSFragger_4_0": "msfragger:hyperscore"},
            # "bigger_scores_better": {"MSFragger_4_0": True},
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


def test_engine_parsers_pglyco_metadata():
    input_file = pytest._test_path / "data" / "test_pglyco_3.txt"
    parser = PGlyco_3_Parser(
        input_file,
        params={
            "cpus": 2,
            "enzyme": "(?<=[KR])(?![P])",
            "terminal_cleavage_site_integrity": "any",
            # "validation_score_field": {"MSFragger_4_0": "msfragger:hyperscore"},
            # "bigger_scores_better": {"MSFragger_4_0": True},
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
    assert parser.metadata["File Origin"] == "pGlyco"


def test_engine_parsers_pglyco_check_parser_compatibility():
    input_file = pytest._test_path / "data" / "test_pglyco_3.txt"
    assert PGlyco_3_Parser.check_parser_compatibility(input_file) is True


def test_engine_parsers_pglyco_check_dataframe_integrity():
    input_file = pytest._test_path / "data" / "test_pglyco_3.txt"
    rt_lookup_path = pytest._test_path / "data" / "test_pglyco_3_meta_data.csv"
    db_path = pytest._test_path / "data" / "test_pglyco_3.fasta"

    parser = PGlyco_3_Parser(
        input_file,
        params={
            "cpus": 2,
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "enzyme": "(?<=[KR])(?![P])",
            "terminal_cleavage_site_integrity": "any",
            # "validation_score_field": {"MSFragger_4_0": "msfragger:hyperscore"},
            # "bigger_scores_better": {"MSFragger_4_0": True},
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
    assert len(df) == 3417
    assert pytest.approx(df["ucalc_mz"].mean()) == 477.8585
    assert pytest.approx(df["exp_mz"].mean()) == 478.12137

    assert df["modifications"].str.contains("Acetyl:0").sum() == 2
    assert df["modifications"].str.contains("Oxidation:").sum() == 221
    assert (
        df["modifications"].str.count("Carbamidomethyl:")
        == df["sequence"].str.count("C")
    ).all()
    assert df["modifications"].str.count(":").sum() == 2242
    assert (df["raw_data_location"] == "path/for/glory.mzML").all()


# def test_pglyco_convert_glycan_composition():
#     input_file = pytest._test_path / "data" / "BSA1_open_search.msfragger4.tsv"
#     rt_lookup_path = pytest._test_path / "data" / "BSA1_ursgal_lookup.csv"
#     db_path = pytest._test_path / "data" / "BSA.fasta"

#     parser = PGlyco_3_Parser(
#         input_file,
#         params={
#             "cpus": 2,
#             "rt_pickle_name": rt_lookup_path,
#             "database": db_path,
#             "enzyme": "(?<=[KR])(?![P])",
#             "terminal_cleavage_site_integrity": "any",
#             # "validation_score_field": {"MSFragger_4_0": "msfragger:hyperscore"},
#             # "bigger_scores_better": {"MSFragger_4_0": True},
#             "modifications": [],
#         },
#     )
#     df = parser.unify()
#     assert df["glycan_composition"].str.contains("Acetyl:0").sum() == 2
