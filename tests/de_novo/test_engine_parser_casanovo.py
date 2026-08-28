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
    assert len(df) == 136
    assert pytest.approx(df["ucalc_mz"].mean(), abs=1e-3) == 455.669
    assert pytest.approx(df["exp_mz"].mean(), abs=1e-3) == 460.018

    assert df["modifications"].str.contains("Carbamyl:0").sum() == 27
    assert df["modifications"].str.contains("Oxidation:").sum() == 29
    assert (
        df["modifications"].str.count("Carbamidomethyl:")
        == df["sequence"].str.count("C")
    ).all()
    assert df["modifications"].str.count(":").sum() == 105
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
    # Manual injection of a dataframe to test specific mapping logic.
    # NOTE: the new parser reads modifications from a separate mzTab-style
    # "modifications" column (e.g. "2-Oxidation (M):UNIMOD:35"), not from
    # bracket notation embedded directly in the sequence string.
    parser.df = pd.DataFrame(
        {
            "sequence": ["LMDKPEQLR"],
            "modifications": ["2-Oxidation (M):UNIMOD:35"],
            "spectrum_id": ["ms_run[1]:controllerType=0 controllerNumber=1 scan=4303"],
            "charge": [2],
            "casanovo:search_engine_score[1]": [0.99],
            "protein_id": [None],
            "exp_mz": [447.34625],
            "retention_time_seconds": [213.85],
        }
    )
    df = parser.unify()
    # Resolved name should carry the same position the mass-shift entry specified.
    assert df["modifications"].iloc[0] == "Oxidation:2"
    assert df["sequence"].iloc[0] == "LMDKPEQLR"


def test_map_mod_translation_casanovo_mass_shift():
    """Mass-only modification (no UNIMOD name in mzTab) should resolve via
    mod_mapper.mass_to_names(), leaving nothing but the mass delta if unresolved.
    """
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

    parser.df = pd.DataFrame(
        {
            "sequence": ["LMDKPEQLR"],
            "modifications": ["2-[15.9949]"],
            "spectrum_id": ["ms_run[1]:controllerType=0 controllerNumber=1 scan=4266"],
            "charge": [2],
            "casanovo:search_engine_score[1]": [0.98],
            "protein_id": [None],
            "exp_mz": [447.34625],
            "retention_time_seconds": [214.10],
        }
    )

    df = parser.unify()

    assert df["modifications"].iloc[0] == "Oxidation:2"
    assert df["sequence"].iloc[0] == "LMDKPEQLR"


def test_map_mod_translation_casanovo_mass_shift_no_mass_to_name(capsys):
    """Mass-only modification (no UNIMOD name in mzTab) should resolve via
    mod_mapper.mass_to_names(), leaving nothing but the mass delta if unresolved.
    """
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

    parser.df = pd.DataFrame(
        {
            "sequence": ["LMDKPEQLR"],
            "modifications": ["0-[+25.980265]"],
            "spectrum_id": ["ms_run[1]:controllerType=0 controllerNumber=1 scan=4266"],
            "charge": [2],
            "casanovo:search_engine_score[1]": [0.98],
            "protein_id": [None],
            "exp_mz": [564.7929],
            "retention_time_seconds": [214.10],
        }
    )
    capsys.readouterr()
    df = parser.unify()

    # Check that warning that mass to names was unsuccessful is dusplayed.
    captured = capsys.readouterr()
    assert "Warning" in captured.out
    assert "25.980265" in captured.out
    assert "no UNIMOD match found" in captured.out
    # Modifications column should be empty.
    assert not df["modifications"].iloc[0] or pd.isna(df["modifications"].iloc[0])
    # Verify that mass_delta is strictly between 25 and 28
    assert 25 < df["mass_delta"].iloc[0] < 28

    assert df["sequence"].iloc[0] == "LMDKPEQLR"


def test_map_mod_translation_casanovo_multiple_mass_shifts():
    """Multiple mass-shift modifications (no UNIMOD names) that cannot be resolved
    should result in an empty modifications column.
    """
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

    parser.df = pd.DataFrame(
        {
            "sequence": ["LMDKPEQLR"],
            "modifications": ["0-[+25.980265];3-[+18.010565]"],
            "spectrum_id": ["ms_run[1]:controllerType=0 controllerNumber=1 scan=4266"],
            "charge": [2],
            "casanovo:search_engine_score[1]": [0.98],
            "protein_id": [None],
            "exp_mz": [578.295],
            "retention_time_seconds": [214.10],
        }
    )
    df = parser.unify()
    # When multiple mass shifts cannot be resolved, modifications column should be empty
    assert not df["modifications"].iloc[0] or pd.isna(df["modifications"].iloc[0])

    assert df["sequence"].iloc[0] == "LMDKPEQLR"
