import pytest

from pyiohat.parsers.ident.ptmshepherd_parser import PTMShepherd_Parser


def test_engine_parsers_ptmshepherd_init():
    input_file = (
        pytest._test_path / "data" / "ptmshepherd_parser_input_file.tsv"
    )
    parser = PTMShepherd_Parser(
        input_file,
        params={
            "cpus": 2,
            "enzyme": "(?<=[KR])(?![P])",
            "terminal_cleavage_site_integrity": "any",
            "validation_score_field": {"MSFragger_4_0": "msfragger:hyperscore"},
            "bigger_scores_better": {"MSFragger_4_0": True},
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


def test_engine_parsers_ptmshepherd_metadata():
    input_file = (
        pytest._test_path / "data" / "ptmshepherd_parser_input_file.tsv"
    )
    parser = PTMShepherd_Parser(
        input_file,
        params={
            "cpus": 2,
            "enzyme": "(?<=[KR])(?![P])",
            "terminal_cleavage_site_integrity": "any",
            "validation_score_field": {"MSFragger_4_0": "msfragger:hyperscore"},
            "bigger_scores_better": {"MSFragger_4_0": True},
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
    assert parser.metadata["validation_score_field"] == "msfragger:hyperscore"
    assert parser.metadata["bigger_scores_better"] == True


def test_engine_parsers_ptmshepherd_check_parser_compatibility():
    input_file = (
        pytest._test_path / "data" / "ptmshepherd_parser_input_file.tsv"
    )
    assert PTMShepherd_Parser.check_parser_compatibility(input_file) is True


def test_engine_parsers_ptmshepherd_check_dataframe_integrity():
    input_file = pytest._test_path / "data" / "BSA1_msfragger_4.tsv"
    rt_lookup_path = pytest._test_path / "data" / "BSA1_ursgal_lookup.csv"
    db_path = pytest._test_path / "data" / "BSA.fasta"

    parser = PTMShepherd_Parser(
        input_file,
        params={
            "cpus": 2,
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "enzyme": "(?<=[KR])(?![P])",
            "terminal_cleavage_site_integrity": "any",
            "validation_score_field": {"MSFragger_4_0": "msfragger:hyperscore"},
            "bigger_scores_better": {"MSFragger_4_0": True},
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