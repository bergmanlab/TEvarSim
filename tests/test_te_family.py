"""
Tests for the family name Simulate writes into INFO/TEFAMILY.

An insertion's ID ends in "_" and its modification string ("TY1-FULL#LTR/Copia_",
"LINE_1SNP0INDEL5polyA"), and everything before that underscore is the family. A deletion's or
an excision's ID ends in the family name itself, which may contain an underscore: splitting
every ID at its last underscore cut "DEL-...-LTR/Gypsy-TY3_1p-FULL" to "...TY3".

No external test runner is required::

    PYTHONPATH=. python tests/test_te_family.py

The test functions are also discoverable by pytest.
"""
import os
import sys

# Import the working-tree package regardless of any installed copy.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from TEvarSim import simulate  # noqa: E402

parse = simulate.Simulator._parse_te_modification


def test_insertion_without_modifications_drops_the_trailing_underscore():
    assert parse(simulate.Simulator, "TY1-FULL#LTR/Copia_") == ("TY1-FULL#LTR/Copia", {})


def test_insertion_family_with_an_underscore_keeps_it():
    assert parse(simulate.Simulator, "SVA_A#Retroposon/SVA_") == ("SVA_A#Retroposon/SVA", {})


def test_insertion_modifications_are_read():
    assert parse(simulate.Simulator, "LINE_1SNP0INDEL5polyA") == ("LINE", {"nSNP": 1, "npolyA": 5})


def test_deletion_family_with_an_underscore_is_not_cut():
    te_id = "DEL-chrIII-96097-101507-LTR/Gypsy-TY3_1p-FULL"
    assert parse(simulate.Simulator, te_id) == (te_id, {})


def test_excision_family_with_an_underscore_is_not_cut():
    te_id = "EXC-chrXIV-611713-617630-335-LTR/Gypsy-TY3_1p-FULL"
    assert parse(simulate.Simulator, te_id) == (te_id, {})


def test_deletion_without_an_underscore_is_unchanged():
    te_id = "DEL-chr21-39156596-39158198-Retroposon/SVA-SVA"
    assert parse(simulate.Simulator, te_id) == (te_id, {})


if __name__ == "__main__":
    failed = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"ok   {name}")
            except AssertionError as e:
                failed += 1
                print(f"FAIL {name}: {e}")
    sys.exit(1 if failed else 0)
