import pytest

from pyiohat.parsers.ident.ptmshepherd_parser import PTMShepherd_Parser


def test_engine_parsers_ptmshepherd_init():
    input_file = pytest._test_path / "data" / "ptmshepherd_parser_input_file.tsv"
    parser = PTMShepherd_Parser(
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


def test_engine_parsers_ptmshepherd_metadata():
    input_file = pytest._test_path / "data" / "ptmshepherd_parser_input_file.tsv"
    parser = PTMShepherd_Parser(
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
    assert parser.metadata["validation_score_field"] == "msfragger:hyperscore"
    assert parser.metadata["bigger_scores_better"] == True


def test_engine_parsers_ptmshepherd_check_parser_compatibility():
    input_file = pytest._test_path / "data" / "ptmshepherd_parser_input_file.tsv"
    assert PTMShepherd_Parser.check_parser_compatibility(input_file) is True
