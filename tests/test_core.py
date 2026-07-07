import json
from pathlib import Path

import pytest

from histo_altloc_cleaner import AltlocSplitter, StructureError, split_altlocs
from histo_altloc_cleaner.core import load_structure

FIXTURE = Path(__file__).parent / "fixtures" / "7s79_1_aligned.cif"
NO_ALTLOC_FIXTURE = Path(__file__).parent / "fixtures" / "1hhk_1_peptide.cif"


def test_altloc_labels():
    splitter = AltlocSplitter(FIXTURE)
    assert splitter.altloc_labels() == ["A", "B"]


def test_summary_lists_all_disordered_residues():
    splitter = AltlocSplitter(FIXTURE)
    summary = splitter.summary()
    assert summary["source"] == FIXTURE.name
    assert summary["altloc_labels"] == ["A", "B"]
    assert len(summary["residues"]) == 23

    residue_131 = next(r for r in summary["residues"] if r["chain"] == "A" and r["residue"] == 131)
    assert residue_131["resname"] == "ARG"
    assert residue_131["disordered_atoms"]["CA"] == {"A": 0.5, "B": 0.5}


def test_split_produces_one_file_per_label_and_a_json_summary(tmp_path):
    splitter = AltlocSplitter(FIXTURE)
    result = splitter.split(tmp_path, fmt="cif")

    assert set(result["labels"]) == {"A", "B"}
    for path in result["labels"].values():
        assert path.exists()
    assert result["summary_path"].exists()

    written = json.loads(result["summary_path"].read_text())
    assert written["outputs"]["A"]["file"] == result["labels"]["A"].name
    assert written["outputs"]["B"]["file"] == result["labels"]["B"].name


def test_split_output_has_no_remaining_disorder(tmp_path):
    splitter = AltlocSplitter(FIXTURE)
    result = splitter.split(tmp_path, fmt="cif")

    for path in result["labels"].values():
        structure = load_structure(path)
        for atom in structure[0].get_atoms():
            assert not atom.is_disordered()


def test_split_selects_correct_conformer_coordinates(tmp_path):
    splitter = AltlocSplitter(FIXTURE)
    result = splitter.split(tmp_path, fmt="cif")

    original = load_structure(FIXTURE)
    residue = original[0]["A"][4]
    original_n = residue["N"]
    coord_a = original_n.child_dict["A"].coord
    coord_b = original_n.child_dict["B"].coord

    out_a = load_structure(result["labels"]["A"])
    out_b = load_structure(result["labels"]["B"])
    assert list(out_a[0]["A"][4]["N"].coord) == pytest.approx(list(coord_a))
    assert list(out_b[0]["A"][4]["N"].coord) == pytest.approx(list(coord_b))


def test_split_blanks_altloc_and_preserves_occupancy(tmp_path):
    splitter = AltlocSplitter(FIXTURE)
    result = splitter.split(tmp_path, fmt="cif")

    original = load_structure(FIXTURE)
    original_occ_a = original[0]["A"][4]["N"].child_dict["A"].get_occupancy()

    out_a = load_structure(result["labels"]["A"])
    atom = out_a[0]["A"][4]["N"]
    assert atom.get_altloc() == " "
    assert atom.get_occupancy() == original_occ_a


def test_split_non_disordered_atoms_identical_across_labels(tmp_path):
    splitter = AltlocSplitter(FIXTURE)
    result = splitter.split(tmp_path, fmt="cif")

    out_a = load_structure(result["labels"]["A"])
    out_b = load_structure(result["labels"]["B"])

    # Residue 1 (chain A) has no altlocs in this fixture; its atoms must
    # be identical in both outputs.
    res_a = out_a[0]["A"][1]
    res_b = out_b[0]["A"][1]
    for atom_name in res_a.child_dict:
        assert list(res_a[atom_name].coord) == pytest.approx(list(res_b[atom_name].coord))
        assert res_a[atom_name].get_occupancy() == res_b[atom_name].get_occupancy()


def test_split_pdb_format_is_parseable(tmp_path):
    splitter = AltlocSplitter(FIXTURE)
    result = splitter.split(tmp_path, fmt="pdb")

    for label, path in result["labels"].items():
        assert path.suffix == ".pdb"
        structure = load_structure(path)
        assert len(list(structure[0].get_atoms())) > 0


def test_no_altlocs_raises():
    with pytest.raises(StructureError):
        AltlocSplitter(NO_ALTLOC_FIXTURE)


def test_module_level_convenience_function(tmp_path):
    result = split_altlocs(FIXTURE, tmp_path, fmt="cif")
    assert set(result["labels"]) == {"A", "B"}
