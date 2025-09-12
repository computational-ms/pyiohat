import pytest

from pyiohat.parsers.ident.peptide_forest_3_parser import peptide_forest_3_Parser


def test_engine_parsers_peptide_forest_init():
    input_file = pytest._test_path / "data" / "test_peptide_forest_3.txt"
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


def test_engine_parsers_peptide_forest_metadata():
    input_file = pytest._test_path / "data" / "test_peptide_forest_3.txt"
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
    assert parser.metadata["File Origin"] == "peptide_forest_3"


def test_engine_parsers_peptide_forest_check_parser_compatibility():
    input_file = pytest._test_path / "data" / "test_peptide_forest_3.txt"
    assert peptide_forest_3_Parser.check_parser_compatibility(input_file) is True


def test_engine_parsers_peptide_forest_check_dataframe_integrity():
    input_file = pytest._test_path / "data" / "test_peptide_forest_3.txt"
    rt_lookup_path = pytest._test_path / "data" / "test_peptide_forest_3_meta_data.csv"
    db_path = pytest._test_path / "data" / "test_glyco_decipher_1.fasta"

    parser = peptide_forest_3_Parser(
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
    assert len(df) == 9
    assert pytest.approx(df["ucalc_mz"].mean()) == 226.733215
    assert pytest.approx(df["exp_mz"].mean()) == 869.482910

    assert df["modifications"].str.contains("Carbamidomethyl:1").sum() == 1
    assert (df["glycan_is_decoy"] == True).sum() == 4
    assert (df["peptide_is_decoy"] == True).sum() == 1
    assert (df["is_decoy"] == True).sum() == 5
    assert df["glycan_composition"].str.contains("NulNAcA").sum() == 8
    assert df["glycan_composition"].str.contains("dHex").sum() == 1
    assert df["protein_id"].str.contains("decoy_").sum() == 1
