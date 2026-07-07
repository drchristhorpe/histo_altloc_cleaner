# CLAUDE.md

Guidance for Claude Code when working in this repository.

## What this is

`histo_altloc_cleaner` splits a structure file (PDB/mmCIF) containing
ALTLOC (alternate conformation) records into one complete structure file
per distinct altloc label, plus a JSON summary describing the altlocs
found. It's a Python library, a Click-based CLI
(`histo-altloc-cleaner`) with Rich console output, and a Claude skill
wrapping that CLI (`skills/histo-altloc-cleaner/`). See
[README.md](README.md) for user-facing usage and
[docs/PLAN.md](docs/PLAN.md) for the design rationale.

## Environment

- Python 3.14, managed with `uv`. Use `uv sync`, `uv run <cmd>`, `uv run pytest`.
- Don't invoke a bare `python`/`pip` — always go through `uv run` /
  `uv add` so the lockfile stays authoritative.

## Layout

```
src/histo_altloc_cleaner/
  core.py   # load_structure(), AltlocSplitter, split_altlocs() convenience wrapper
  cli.py    # Click CLI + Rich table output (entry point: histo-altloc-cleaner)
tests/
  fixtures/  # real structure files used by the tests, keep committed
  test_core.py
  test_cli.py
skills/histo-altloc-cleaner/SKILL.md
```

## Key invariants — don't break these

- `AltlocSplitter(path)` parses the structure **once** in `__init__` for
  *analysis* (detecting altloc labels, occupancies, building the
  fallback choice maps) — but `split()` **re-parses the source file
  fresh for every label's write**. This is deliberate, not an
  inefficiency to "optimize away": `Bio.PDB`'s `get_unpacked_list()`
  returns the actual child `Atom` objects held by a `DisorderedAtom`
  container (not copies), and the write path mutates the accepted atom's
  `.altloc` in place (see below). Sharing one parsed `Structure` across
  multiple labels' writes would let label `A`'s write corrupt data label
  `B`'s write still needs.
- Only atom-level disorder (`DisorderedAtom`) is handled. Whole-residue
  microheterogeneity (`DisorderedResidue` — different residue *types* at
  the same position) is out of scope; it isn't exercised by the
  committed fixture and hasn't come up. Don't assume it's handled without
  checking first.
- **Mismatched altloc labels**: if a disordered position doesn't carry
  the requested output label (e.g. the structure's global label set is
  `A`/`B` but one residue is only disordered as `A`/`C`), the atom falls
  back to whichever conformer at that position has the **highest
  occupancy** (ties broken alphabetically by code). This is an explicit
  choice, confirmed with the user, so every output file stays structurally
  complete. Every fallback is recorded in the JSON summary's
  `outputs.<label>.fallback_atoms` — never silent.
- **Output altloc/occupancy handling**: written atoms at a previously
  disordered position have their altloc code **blanked** (`" "` in
  PDB / `.` in mmCIF), since the file now represents one complete
  conformer, not a mixture. **Occupancy is left untouched** — the
  original value is preserved as data provenance, not reset to `1.0`.
  This was an explicit choice, confirmed with the user; don't change
  either behavior without checking first.
- The JSON summary describes the **source** structure's altlocs (one
  summary per run, not one per output file) — `residues` lists every
  disordered residue found, `outputs` maps each label to its output
  filename and fallback list.
- Only the first model (`structure[0]`) is used everywhere — sufficient
  for X-ray/cryo-EM, a documented limitation for NMR ensembles.
- A structure with zero disordered atoms raises `StructureError` — there's
  nothing to split.

## Testing

- `uv run pytest` — fixtures are already committed in `tests/fixtures/`;
  tests don't hit the network.
- `7s79_1_aligned.cif` is the primary fixture: a real deposited structure
  (downloaded from `coordinates.histo.fyi`) with genuine ALTLOC `A`/`B`
  records — 23 disordered residues across two chains, varied occupancy
  splits (0.5/0.5, 0.7/0.3, 0.8/0.2). Unlike the `8gvi_1_aligned` fixture
  used by sibling tools, whose alignment pipeline strips altlocs
  entirely, this entry retains them — verified by inspection before
  choosing it.
- `1hhk_1_peptide.cif` (copied from `histo_com`'s fixtures) has no
  altlocs at all and is used only for the `StructureError` "no ALTLOCs
  found" error-path test — reused rather than hand-writing a synthetic
  structure, per the project-wide real-fixtures-only rule.
- The mismatched-altloc-label fallback path (§ above) isn't covered by a
  dedicated test — no committed real fixture exercises it, and per the
  real-fixtures-only convention a synthetic structure wasn't created just
  to hit that branch. If a fixture surfaces this case, add a test for it.

## Scope

The CLI intentionally exposes exactly two options: `--output-dir` and
`--format` (plus the positional `FILENAME`). Don't add further options
(e.g. explicit label selection, occupancy-reset toggle, whole-residue
disorder handling) without checking with the user first — these were
deliberate constraints from the initial design conversation, not
oversights.
