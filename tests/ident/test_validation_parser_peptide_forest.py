import pytest

from pyiohat.parsers.ident.peptide_forest_parser import PeptideForest_Parser


def test_engine_parsers_peptide_forest_init():
    input_file = pytest._test_path / "data" / "test_peptide_forest.csv"
    parser = PeptideForest_Parser(
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
                    "name": "Methylthio",
                },
                {
                    "aa": "K",
                    "type": "opt",
                    "position": "any",
                    "name": "Label:13C(6)15N(2)",
                },
            ],
        },
    )


def test_engine_parsers_peptide_forest_metadata():
    input_file = pytest._test_path / "data" / "test_peptide_forest.csv"
    parser = PeptideForest_Parser(
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
                    "name": "Methylthio",
                },
                {
                    "aa": "K",
                    "type": "opt",
                    "position": "any",
                    "name": "Label:13C(6)15N(2)",
                },
            ],
        },
    )
    assert parser.metadata["File Origin"] == "peptide_forest_3"


def test_engine_parsers_peptide_forest_check_parser_compatibility():
    input_file = pytest._test_path / "data" / "test_peptide_forest.csv"
    assert PeptideForest_Parser.check_parser_compatibility(input_file) is True


def test_engine_parsers_peptide_forest_check_dataframe_integrity():
    input_file = pytest._test_path / "data" / "test_peptide_forest.csv"
    rt_lookup_path = pytest._test_path / "data" / "test_peptide_forest_meta_data.csv"
    db_path = pytest._test_path / "data" / "Hfvol_prot_250410.fasta"

    parser = PeptideForest_Parser(
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
                    "name": "Methylthio",
                },
                {
                    "aa": "K",
                    "type": "opt",
                    "position": "any",
                    "name": "Label:13C(6)15N(2)",
                },
            ],
        },
    )
    df = parser.unify()
    assert len(df) == 9
    assert pytest.approx(df["ucalc_mz"].mean()) == 435.864868
    assert pytest.approx(df["exp_mz"].mean()) == 435.8951
    print(df["sequence_stop"])
    # assert df["modifications"].str.contains("Label:13C(6)15N(2):12").sum() == 5
    assert df["sequence_start"].str.contains("23").sum() == 8
    assert df["sequence_stop"].str.contains("475").sum() == 1
