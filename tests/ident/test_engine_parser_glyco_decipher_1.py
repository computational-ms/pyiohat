import pytest

from pyiohat.parsers.ident.glyco_decipher_1_parser import GlycoDecipher_1_Parser


def test_engine_parsers_glyco_decipher_init():
    input_file = pytest._test_path / "data" / "test_glyco_decipher_1.txt"
    parser = GlycoDecipher_1_Parser(
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


def test_engine_parsers_glyco_decipher_metadata():
    input_file = pytest._test_path / "data" / "test_glyco_decipher_1.txt"
    parser = GlycoDecipher_1_Parser(
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
    assert parser.metadata["File Origin"] == "GlycoDecipher"


def test_engine_parsers_glyco_decipher_check_parser_compatibility():
    input_file = pytest._test_path / "data" / "test_glyco_decipher_1.txt"
    assert GlycoDecipher_1_Parser.check_parser_compatibility(input_file) is True


def test_engine_parsers_glyco_decipher_check_dataframe_integrity():
    input_file = pytest._test_path / "data" / "test_glyco_decipher_1.txt"
    rt_lookup_path = pytest._test_path / "data" / "test_glyco_decipher_1_meta_data.csv"
    db_path = pytest._test_path / "data" / "test_glyco_decipher_1.fasta"

    parser = GlycoDecipher_1_Parser(
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
    assert len(df) == 3
    assert pytest.approx(df["ucalc_mz"].mean()) == 322.7364
    assert pytest.approx(df["exp_mz"].mean()) == 1136.3010

    assert df["modifications"].str.contains("Carbamidomethyl:4").sum() == 1
    assert (df["is_decoy"] == True).sum() == 0
    assert df["glycan_composition"].str.contains("NulNAcA").sum() == 3
    assert df["glycan_composition"].str.contains("dHex").sum() == 3
