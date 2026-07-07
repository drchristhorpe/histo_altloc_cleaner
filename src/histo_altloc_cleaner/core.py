"""Structure loading and ALTLOC splitting."""

from __future__ import annotations

import json
from pathlib import Path

from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB.PDBIO import PDBIO, Select
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Structure import Structure

_CIF_SUFFIXES = {".cif", ".mmcif"}
_PDB_SUFFIXES = {".pdb", ".ent"}

_WRITERS = {"pdb": PDBIO, "cif": MMCIFIO}


class StructureError(ValueError):
    """Raised for problems loading or splitting a structure."""


def load_structure(path: str | Path, structure_id: str | None = None) -> Structure:
    """Parse a PDB or mmCIF file into a Bio.PDB Structure.

    Format is chosen from the file extension (case-insensitive):
    ``.cif``/``.mmcif`` -> mmCIF, ``.pdb``/``.ent`` -> legacy PDB.
    """
    path = Path(path)
    if not path.is_file():
        raise StructureError(f"No such file: {path}")

    suffix = path.suffix.lower()
    sid = structure_id or path.stem

    if suffix in _CIF_SUFFIXES:
        parser = MMCIFParser(QUIET=True)
    elif suffix in _PDB_SUFFIXES:
        parser = PDBParser(QUIET=True)
    else:
        raise StructureError(
            f"Unrecognised structure file extension {suffix!r} for {path}; "
            "expected one of .cif, .mmcif, .pdb, .ent"
        )

    structure = parser.get_structure(sid, str(path))
    if len(structure) == 0:
        raise StructureError(f"No models found in {path}")
    return structure


def _disorder_map(
    model,
) -> tuple[dict[tuple[str, tuple, str], dict[str, float]], dict[tuple[str, tuple], str]]:
    """Map ``(chain_id, residue_id, atom_name) -> {altloc_code: occupancy}``
    for every disordered atom position in ``model``, plus a
    ``(chain_id, residue_id) -> resname`` map for the residues involved."""
    positions: dict[tuple[str, tuple, str], dict[str, float]] = {}
    resnames: dict[tuple[str, tuple], str] = {}
    for chain in model:
        for residue in chain:
            for atom in residue:
                if not atom.is_disordered():
                    continue
                key = (chain.id, residue.id, atom.get_name())
                positions[key] = {
                    code: child.get_occupancy()
                    for code, child in atom.child_dict.items()
                }
                resnames[(chain.id, residue.id)] = residue.resname
    return positions, resnames


class AltlocSplitter:
    """Loads a structure once and splits it into one complete structure
    file per distinct ALTLOC label found.

    >>> splitter = AltlocSplitter("structure.cif")
    >>> splitter.altloc_labels()
    ['A', 'B']
    >>> splitter.split("out/", fmt="cif")
    """

    def __init__(self, path: str | Path):
        self.path = Path(path)
        structure = load_structure(self.path)
        self._positions, self._resnames = _disorder_map(structure[0])
        if not self._positions:
            raise StructureError(f"No ALTLOC records found in {self.path}")

    def altloc_labels(self) -> list[str]:
        """Every distinct altloc code found in the structure, sorted."""
        labels: set[str] = set()
        for occupancies in self._positions.values():
            labels.update(occupancies)
        return sorted(labels)

    def summary(self) -> dict:
        """A JSON-serialisable description of the ALTLOCs found, grouped
        by residue (``chain``, ``residue``, ``insertion_code``, ``resname``,
        ``disordered_atoms``: atom name -> {altloc_code: occupancy})."""
        residues: dict[tuple[str, tuple], dict] = {}
        for (chain_id, res_id, atom_name), occupancies in self._positions.items():
            _, resseq, icode = res_id
            res_key = (chain_id, res_id)
            if res_key not in residues:
                residues[res_key] = {
                    "chain": chain_id,
                    "residue": resseq,
                    "insertion_code": icode.strip(),
                    "resname": self._resnames[res_key],
                    "disordered_atoms": {},
                }
            residues[res_key]["disordered_atoms"][atom_name] = dict(occupancies)

        ordered = sorted(residues.values(), key=lambda r: (r["chain"], r["residue"]))
        return {
            "source": self.path.name,
            "altloc_labels": self.altloc_labels(),
            "residues": ordered,
        }

    def _choice_map(self, label: str) -> tuple[dict[tuple, str], list[dict]]:
        """For output ``label``, resolve every disordered position to a
        single altloc code: ``label`` itself if present there, otherwise
        the highest-occupancy conformer (ties broken alphabetically).
        Returns the choice map plus a list of fallback events."""
        choices: dict[tuple, str] = {}
        fallbacks: list[dict] = []
        for key, occupancies in self._positions.items():
            if label in occupancies:
                choices[key] = label
                continue
            fallback_code = max(sorted(occupancies), key=lambda code: occupancies[code])
            choices[key] = fallback_code
            chain_id, res_id, atom_name = key
            fallbacks.append(
                {
                    "chain": chain_id,
                    "residue": res_id[1],
                    "atom": atom_name,
                    "used_altloc": fallback_code,
                    "reason": f"no {label} conformer present",
                }
            )
        return choices, fallbacks

    def split(self, output_dir: str | Path, fmt: str = "cif") -> dict:
        """Write one structure file per altloc label into ``output_dir``
        (created if missing), plus one JSON summary. ``fmt`` is ``"pdb"``
        or ``"cif"``.

        Returns ``{"labels": {label: Path}, "summary_path": Path,
        "summary": dict}``.
        """
        if fmt not in _WRITERS:
            raise ValueError(f"fmt must be one of {sorted(_WRITERS)}, got {fmt!r}")

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        stem = self.path.stem

        summary = self.summary()
        outputs: dict[str, dict] = {}
        label_paths: dict[str, Path] = {}

        for label in self.altloc_labels():
            choices, fallbacks = self._choice_map(label)
            out_path = output_dir / f"{stem}_altloc{label}.{fmt}"
            _write_label(self.path, label, choices, fmt, out_path)
            label_paths[label] = out_path
            outputs[label] = {"file": out_path.name, "fallback_atoms": fallbacks}

        summary["outputs"] = outputs
        summary_path = output_dir / f"{stem}_altlocs.json"
        summary_path.write_text(json.dumps(summary, indent=2))

        return {"labels": label_paths, "summary_path": summary_path, "summary": summary}


class _AltlocSelect(Select):
    """Keeps every ordered atom, and exactly one conformer per disordered
    position (per ``choices``), blanking its altloc code."""

    def __init__(self, choices: dict[tuple, str], label: str):
        self.choices = choices
        self.label = label

    def accept_atom(self, atom) -> bool:
        if not atom.is_disordered():
            return True
        residue = atom.get_parent()
        chain = residue.get_parent()
        key = (chain.id, residue.id, atom.get_name())
        target = self.choices.get(key, self.label)
        if atom.get_altloc() != target:
            return False
        atom.altloc = " "
        return True


def _write_label(
    source_path: Path, label: str, choices: dict[tuple, str], fmt: str, out_path: Path
) -> None:
    """Re-parse ``source_path`` fresh and write the ``label`` conformer.

    A fresh re-parse (rather than reusing/copying an in-memory Structure)
    is deliberate: ``get_unpacked_list()`` returns the actual child atoms
    of a DisorderedAtom, and ``_AltlocSelect`` mutates the accepted atom's
    altloc in place — reusing one parsed structure across labels would let
    one label's write corrupt the data the next label's write still needs.
    """
    structure = load_structure(source_path)
    io = _WRITERS[fmt]()
    io.set_structure(structure)
    io.save(str(out_path), select=_AltlocSelect(choices, label))


def split_altlocs(path: str | Path, output_dir: str | Path, fmt: str = "cif") -> dict:
    """Convenience wrapper: split a structure file's ALTLOCs into
    ``output_dir`` without holding onto an `AltlocSplitter` instance."""
    return AltlocSplitter(path).split(output_dir, fmt=fmt)
