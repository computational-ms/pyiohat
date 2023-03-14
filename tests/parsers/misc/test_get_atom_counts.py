import numpy as np

from pyiohat.parsers.misc import get_atom_counts


def test_simple():
    compositions = {  # these include water ...
        "E": {"C": 5, "H": 9, "N": 1, "O": 4},
        "Magic": {"H": 100},
    }
    modifications = np.array(
        [
            "",
            "",
            "Magic:1",
        ],
        dtype=str,
    )
    sequences = np.array(
        [
            "E",
            "EE",
            "E",
        ],
        dtype=str,
    )
    elements, atom_counts = get_atom_counts(
        sequences=sequences,
        modifications=modifications,
        compositions=compositions,
    )

    assert np.array_equal(elements, ["C", "H", "N", "O"])
    assert np.array_equal(atom_counts[0], [5, 9, 1, 4])
    assert np.array_equal(atom_counts[1], [10, 16, 2, 7])
    assert np.array_equal(atom_counts[2], [5, 109, 1, 4])


def test_negative():
    compositions = {  # these include water ...
        "E": {"C": 5, "H": 9, "N": 1, "O": 4},
        "Magic": {"H": -10},
    }
    modifications = np.array(
        [
            "Magic:1",
        ],
        dtype=str,
    )
    sequences = np.array(
        [
            "E",
        ],
        dtype=str,
    )
    elements, atom_counts = get_atom_counts(
        sequences=sequences,
        modifications=modifications,
        compositions=compositions,
    )
    assert np.array_equal(atom_counts[0], [5, -1, 1, 4])
