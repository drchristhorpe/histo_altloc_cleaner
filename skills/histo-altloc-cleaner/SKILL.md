---
name: histo-altloc-cleaner
description: Split a structure file (PDB/mmCIF) containing ALTLOC (alternate conformation) records into one complete structure file per distinct altloc, plus a JSON summary of the ALTLOCs found. Use when asked to split, separate, or clean up a structure with alternate conformations/altlocs, or to generate one structure per conformer from a disordered structure file.
---

# histo-altloc-cleaner

`histo-altloc-cleaner` is a CLI tool (installed from the
`histo_altloc_cleaner` package) that splits a PDB/mmCIF structure
containing ALTLOC records into one complete structure file per distinct
altloc label, using Biopython. Invoke it with the Bash tool.

## When to use this skill

The user provides (or references) a `.cif`/`.mmcif` or `.pdb`/`.ent`
structure file that has alternate conformations (ALTLOC records) and asks
to split it into separate files per conformer, "clean up" the altlocs, or
generate one structure file per alternate location.

## Checking availability

```bash
histo-altloc-cleaner --help
```

If this fails with "command not found", install it first:

```bash
uv tool install histo_altloc_cleaner   # or: pip install histo_altloc_cleaner
```

(If working from a checkout of the `histo_altloc_cleaner` source repo
instead of an installed package, use `uv run histo-altloc-cleaner ...`
there instead.)

## Usage

```
histo-altloc-cleaner FILENAME --output-dir DIR [--format pdb|cif]
```

- `FILENAME`: the structure file to split.
- `--output-dir`/`-o` (required): directory to write the split structure
  files and JSON summary into (created if missing).
- `--format`/`-f` (optional, default `cif`): output structure format,
  `pdb` or `cif`. Independent of the input file's format.

Example: `histo-altloc-cleaner structure.cif --output-dir out/`

If the structure has 2 ALTLOCs (e.g. `A`/`B`), this produces
`out/structure_altlocA.cif`, `out/structure_altlocB.cif`, and
`out/structure_altlocs.json`. A structure with no ALTLOC records errors
out — there's nothing to split.

## Output

- One structure file per altloc label (`<stem>_altloc<LABEL>.<ext>`):
  every ordinary atom is present unchanged; every disordered position is
  resolved to its `<LABEL>` conformer. The written altloc code is blanked
  (the file now represents one complete conformer); occupancy values are
  left as-is.
- One JSON summary (`<stem>_altlocs.json`) describing every disordered
  residue found in the source file (chain, residue, resname, and the
  occupancy of each altloc-bearing atom), plus, per output label, the
  filename written and any atoms that had to fall back to a
  different-than-requested conformer (when a position doesn't carry that
  label — recorded, never silently dropped).

The CLI also prints a Rich table (one row per label: output file,
disordered-residue count, fallback-atom count) to the console.

## Example

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

Report the result back to the user in whatever form they asked for (the
list of output files, a residue-count summary, the full JSON, etc.) —
this skill only tells you how to obtain the split structures. Only
atom-level disorder is handled (not whole-residue microheterogeneity),
and only the structure's first model is used.
