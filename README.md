# histo-altloc-cleaner

Split a 3D biological structure file (PDB/mmCIF) that contains **ALTLOC**
(alternate conformation) records into one complete structure file per
distinct altloc — a structure with disordered residues labelled `A`/`B`
produces two output files, each with every disordered atom resolved to a
single conformer, plus a JSON file describing the ALTLOCs found.

Built on [Biopython](https://biopython.org/), it ships as:

- a Python library — `import histo_altloc_cleaner`
- a CLI tool — `histo-altloc-cleaner`
- a [Claude Code / Claude Desktop skill](skills/histo-altloc-cleaner/SKILL.md)

Requires Python 3.14+.

## Install

```bash
uv sync                 # dev environment, from a checkout
uv tool install .       # install the `histo-altloc-cleaner` CLI globally
# or
pip install .
```

## CLI usage

```
histo-altloc-cleaner FILENAME --output-dir DIR [--format pdb|cif]
```

- `FILENAME` — a `.cif`/`.mmcif` or `.pdb`/`.ent` structure file
  containing ALTLOC records.
- `--output-dir`, `-o` — directory to write the split structure files and
  JSON summary into (created if missing).
- `--format`, `-f` — output structure format, `pdb` or `cif` (default
  `cif`). Independent of the input file's format.

### Example

```bash
$ histo-altloc-cleaner 7s79_1_aligned.cif --output-dir out/
                ALTLOC split: 7s79_1_aligned.cif
┏━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━┓
┃ label ┃ output file                   ┃ disordered residues ┃ fallback atoms ┃
┡━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━┩
│ A     │ out/7s79_1_aligned_altlocA.cif │ 23                  │ 0              │
│ B     │ out/7s79_1_aligned_altlocB.cif │ 23                  │ 0              │
└───────┴───────────────────────────────┴─────────────────────┴────────────────┘
2 altloc(s) (A, B) written; summary at out/7s79_1_aligned_altlocs.json
```

This writes `out/7s79_1_aligned_altlocA.cif`,
`out/7s79_1_aligned_altlocB.cif`, and
`out/7s79_1_aligned_altlocs.json`.

## Output files

- One structure file per altloc label: `<stem>_altloc<LABEL>.<ext>`. Every
  non-disordered atom is present, unchanged, in every output file; every
  disordered position is resolved to its `<LABEL>` conformer (or, if that
  position doesn't carry `<LABEL>`, the highest-occupancy conformer
  present there — see "Notes and limitations"). The written altloc code
  is blanked (since the file now represents one complete conformer), but
  the original occupancy value is left untouched.
- One JSON summary, `<stem>_altlocs.json`, describing the ALTLOCs found in
  the *source* file:

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

`fallback_atoms` lists any disordered atom that didn't carry the output
label and had to fall back to another conformer, e.g. `{"chain": "A",
"residue": 61, "atom": "OG", "used_altloc": "A", "reason": "no B
conformer present"}`.

## Library usage

```python
from histo_altloc_cleaner import AltlocSplitter

splitter = AltlocSplitter("7s79_1_aligned.cif")
splitter.altloc_labels()              # -> ["A", "B"]
splitter.summary()                    # -> dict (see JSON shape above, minus "outputs")
result = splitter.split("out/", fmt="cif")
# -> {"labels": {"A": Path("out/7s79_1_aligned_altlocA.cif"), "B": Path(...)},
#     "summary_path": Path("out/7s79_1_aligned_altlocs.json"),
#     "summary": {...}}
```

A one-shot convenience function is also available:

```python
from histo_altloc_cleaner import split_altlocs

split_altlocs("7s79_1_aligned.cif", "out/", fmt="cif")
```

## Notes and limitations

- Only the **first model** in a file is used — sufficient for
  X-ray/cryo-EM structures; NMR ensembles are not averaged/iterated
  across models.
- Only atom-level disorder is handled (the common case: a residue's
  sidechain or backbone atoms split into alternate positions). Whole
  -residue microheterogeneity (different residue *types* modelled at the
  same position) is out of scope.
- If a disordered position doesn't carry the output label being written
  (e.g. the file's global label set is `A`/`B` but one residue is only
  disordered as `A`/`C`), the highest-occupancy conformer present at that
  position is used instead, and the substitution is recorded in the JSON
  summary's `fallback_atoms` list — nothing is silently dropped.
- A structure with no ALTLOC records raises an error — there's nothing to
  split.

## Development

```bash
uv sync
uv run pytest
```

Test fixtures under `tests/fixtures/` are real structure files:
`7s79_1_aligned.cif` (from
[coordinates.histo.fyi](https://coordinates.histo.fyi/), genuine ALTLOC
`A`/`B` records) and `1hhk_1_peptide.cif` (no altlocs, used for the
error-path test).

See [docs/PLAN.md](docs/PLAN.md) for the design rationale and
[CHANGELOG.md](CHANGELOG.md) for release history.
