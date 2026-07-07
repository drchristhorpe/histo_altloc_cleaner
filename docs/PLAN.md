# histo_altloc_cleaner — Design & Implementation Plan

## 1. Purpose

`histo_altloc_cleaner` splits a 3D biological structure file (PDB/mmCIF)
that contains **ALTLOC** (alternate conformation) records into one
complete structure file per distinct altloc label — e.g. a structure with
disordered residues carrying labels `A`/`B` produces two output structure
files, one with every disordered atom resolved to its `A` conformer, one
with every disordered atom resolved to its `B` conformer. Non-disordered
atoms are identical (and present) in every output file. Alongside the
structure files, it writes one JSON file describing the altlocs found in
the source structure. It ships as:

1. A Python library (`import histo_altloc_cleaner`)
2. A CLI tool (`histo-altloc-cleaner`), built with Click, with
   Rich-formatted console output
3. A Claude Code / Claude Desktop skill that wraps the CLI

## 2. Tooling

- Python **3.14**, managed with **uv** (`uv venv`, `uv sync`, `uv run`, `uv build`)
- **Biopython** for structure parsing (`Bio.PDB.MMCIFParser`/`PDBParser`)
  and writing (`Bio.PDB.PDBIO`/`Bio.PDB.MMCIFIO`)
- **Click** for the CLI, **Rich** for console summary output
- `pyproject.toml` with a `src/` layout, `uv_build` build backend, and a
  `[project.scripts]` entry point (`histo-altloc-cleaner = "histo_altloc_cleaner.cli:main"`)
- `pytest` for tests

## 3. Structure loading

`histo_altloc_cleaner.core.load_structure(path)` follows the same pattern
established in `histo_com`/`histo_neighbours`: pick
`MMCIFParser`/`PDBParser` from the file extension (`.cif`/`.mmcif` vs
`.pdb`/`.ent`, case-insensitive), `QUIET=True`, operate only on the
**first model** (`structure[0]`) — sufficient for X-ray/cryo-EM, a
documented limitation for NMR ensembles.

## 4. Altloc detection

`AltlocSplitter(path)` parses the structure once (for analysis only —
see §6 for why writing re-parses) and walks every residue's
`child_dict`: an entry that is a `Bio.PDB.Atom.DisorderedAtom` (rather
than a plain `Atom`) represents a disordered atom position. For each such
position it records, keyed by `(chain_id, residue_id, atom_name)`:

- the set of altloc codes present (`disordered_get_id_list()`)
- each code's occupancy (`child.get_occupancy()`)

The structure's global altloc label set (`altloc_labels()`) is the sorted
union of every code seen at every disordered position, excluding the
blank/ordinary code. This is the set of output files to produce — "2
ALTLOCs found" means 2 output structure files.

A structure with no disordered positions raises `StructureError` ("no
ALTLOC records found") — there's nothing to split.

## 5. Per-label atom selection & the mismatched-label edge case

For a given output label `L` (e.g. `"B"`), every disordered position must
resolve to exactly one atom:

- If `L` is one of the codes present at that position, use it directly.
- If not (a residue disordered as e.g. `A`/`C` when the file's global
  label set is `A`/`B`, so `B` doesn't apply there) — an explicit design
  choice, confirmed with the user: **fall back to the conformer with the
  highest occupancy** at that position (ties broken by alphabetically
  first code, for determinism). This keeps every output file structurally
  complete rather than silently missing atoms. Every fallback is recorded
  in the JSON summary's `outputs.<label>.fallback_atoms` list so it's
  auditable, not silent.

This per-label choice map (`(chain_id, residue_id, atom_name) -> code`)
is computed once from the analysis pass and reused for every label's
write pass.

## 6. Writing output structures

Per the confirmed design: **each label's output is written from a fresh
re-parse of the source file**, not a shared/mutated in-memory structure.
This is deliberate, not an inefficiency to "fix" later: `Bio.PDB`'s
`get_unpacked_list()` returns the *actual* child `Atom` objects held by a
`DisorderedAtom` container (not copies), and blanking the altloc code
post-selection (see below) mutates those objects in place — reusing one
parsed `Structure` across multiple label passes would have label `A`'s
write corrupt the data label `B`'s write still needs. Re-parsing per
label is simple, correct, and cheap at the file sizes this tool targets.

Writing uses `Bio.PDB.PDBIO`/`Bio.PDB.MMCIFIO`'s standard `Select`
mechanism, following the classic Biopython "NotDisordered" altloc-filter
recipe: both writers iterate `residue.get_unpacked_list()` internally, so
`Select.accept_atom(atom)` sees every conformer atom individually with
its own `.get_altloc()`/`.get_occupancy()`.

```python
class _AltlocSelect(Select):
    def accept_atom(self, atom):
        if not atom.is_disordered():
            return True
        key = (chain_id, atom.get_parent().id, atom.get_name())
        if atom.get_altloc() != choice_map.get(key, label):
            return False
        atom.altloc = " "   # blank: this file has one conformer, not an alternate
        return True
```

Per the confirmed design: the written altloc code is **blanked** (`" "`
in PDB, `.` in mmCIF) since each output file now represents a single,
complete conformer rather than a mixture — matching common altloc-split
tools (e.g. pdb-tools' `pdb_selaltloc`). **Occupancy is left untouched** —
the original value is preserved as data provenance, not reset to `1.0`.

## 7. Output format

`--format {pdb,cif}`, default **`cif`** (no 80-column/hybrid-36 limits,
lossless) — independent of the input file's format; a `.pdb` input can be
split into `.cif` outputs and vice versa.

## 8. Output naming & the JSON summary

Given `--output-dir DIR` and an input file with stem `<stem>` (e.g.
`7s79_1_aligned`):

- One structure file per label: `DIR/<stem>_altloc<LABEL>.<ext>`, e.g.
  `7s79_1_aligned_altlocA.cif`, `7s79_1_aligned_altlocB.cif`.
- One JSON summary: `DIR/<stem>_altlocs.json`, describing the ALTLOCs
  found in the *source* structure (not one summary per output file):

```json
{
  "source": "7s79_1_aligned.cif",
  "altloc_labels": ["A", "B"],
  "residues": [
    {
      "chain": "A",
      "residue": 131,
      "insertion_code": "",
      "resname": "ARG",
      "disordered_atoms": {
        "CG": {"A": 0.5, "B": 0.5},
        "CD": {"A": 0.5, "B": 0.5}
      }
    }
  ],
  "outputs": {
    "A": {"file": "7s79_1_aligned_altlocA.cif", "fallback_atoms": []},
    "B": {"file": "7s79_1_aligned_altlocB.cif", "fallback_atoms": []}
  }
}
```

`fallback_atoms` entries (see §5) look like `{"chain": "A", "residue":
61, "atom": "OG", "used_altloc": "A", "reason": "no B conformer present"}`.

## 9. Library API (`histo_altloc_cleaner/core.py`)

```python
from histo_altloc_cleaner import AltlocSplitter

splitter = AltlocSplitter("structure.cif")
splitter.altloc_labels()              # -> ["A", "B"]
splitter.summary()                    # -> dict, JSON-serialisable (see §8, minus "outputs")
result = splitter.split("out/", fmt="cif")
# -> {"labels": {"A": Path("out/structure_altlocA.cif"), "B": Path(...)},
#     "summary_path": Path("out/structure_altlocs.json"),
#     "summary": {...}}   # the full dict written to summary_path, incl. "outputs"
```

Module-level convenience function `split_altlocs(path, output_dir, fmt="cif")`
wraps `AltlocSplitter` for one-shot scripting use.

## 10. CLI (`histo_altloc_cleaner/cli.py`)

```
histo-altloc-cleaner FILENAME --output-dir DIR [--format pdb|cif]
```

- `FILENAME`: positional argument, structure file path
- `--output-dir`/`-o`: required, directory to write output structure
  files + JSON summary into (created if missing)
- `--format`/`-f`: choice `pdb`/`cif`, default `cif`, output structure format

After splitting, the CLI prints a Rich table (one row per altloc label:
label, output path, disordered-residue count, fallback-atom count) plus a
one-line summary, then confirms the JSON summary path.

## 11. Claude skill

`skills/histo-altloc-cleaner/SKILL.md` — describes when/how to invoke the
`histo-altloc-cleaner` CLI (splitting a structure with alternate
conformations into one file per altloc), so Claude Code/Desktop can call
it via the Bash tool once installed.

## 12. Package layout

```
histo_altloc_cleaner/
  pyproject.toml
  README.md
  CLAUDE.md
  CHANGELOG.md
  docs/PLAN.md
  .gitignore
  src/histo_altloc_cleaner/
    __init__.py
    core.py          # load_structure(), AltlocSplitter, split_altlocs()
    cli.py            # Click CLI (entry point: histo-altloc-cleaner)
    py.typed
  skills/histo-altloc-cleaner/SKILL.md
  tests/
    fixtures/
      7s79_1_aligned.cif   # real deposited structure with genuine ALTLOC A/B records
    test_core.py
    test_cli.py
  tmp/
    .gitkeep
```

## 13. Test fixture

`7s79_1_aligned.cif`, downloaded directly from
`https://coordinates.histo.fyi/structures/downloads/class_i/without_solvent/7s79_1_aligned.cif`
per the user's pointer. Unlike the `8gvi_1_aligned` fixture used by
`histo_neighbours`/`histo_com` (whose alignment pipeline strips altlocs
entirely — confirmed by inspection before choosing this file), this entry
retains genuine alternate conformations: 23 disordered residues across
two chains (16 in chain `A`, 7 in chain `B`), altloc codes `A`/`B` only,
with varied occupancy splits (0.5/0.5, 0.7/0.3, 0.8/0.2) — good coverage
for both the common case and occupancy-based fallback logic.

(RCSB's current deposited `1CRN` — the other candidate discussed — was
checked and no longer carries altloc records in its current PDB-format
export, so it wasn't used.)

## 14. Testing plan

Unit tests with pytest against the committed `7s79_1_aligned.cif` fixture:

- `altloc_labels()` returns exactly `["A", "B"]`.
- `summary()` lists all 23 disordered residues with correct per-atom
  occupancies for a sample of known positions (e.g. chain A residue 131).
- `split()` produces exactly 2 structure files + 1 JSON summary; each
  structure file, when re-parsed, has **zero** remaining disordered atoms
  (every position resolved to a single conformer) and the same total atom
  *count* pattern as the other label's output (same positions covered).
- For a spot-checked disordered atom, the `A`-output file's coordinates
  match the original `A`-conformer coordinates (and `B`-output matches
  `B`), not swapped.
- Written atoms at previously-disordered positions have a blank altloc
  code in the output; occupancy is unchanged from the source value.
- Non-disordered atoms are byte-identical (coordinates/occupancy) across
  both output files.
- `--format cif` vs `--format pdb` both produce parseable output for the
  same split.
- Error path: a structure with no altlocs raises `StructureError`.
- CLI test (`click.testing.CliRunner`): writes the expected number of
  files to `--output-dir`, non-zero exit when the file is missing.

## 15. Workflow

1. Write this plan.
2. Scaffold project with `uv init`/`uv add`.
3. Implement `core.py`, `cli.py`.
4. Write `CHANGELOG.md` as work proceeds.
5. Write and run tests against the committed fixture.
6. Write `README.md`, `CLAUDE.md`, the skill file.
7. Manually exercise the CLI end-to-end (`uv run histo-altloc-cleaner
   ...`), inspect the output structure files + JSON in `tmp/`.
8. Present for approval, then commit.
