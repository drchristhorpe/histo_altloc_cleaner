# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

## [0.1.0] - 2026-07-07

### Added

- `histo_altloc_cleaner` Python library: `AltlocSplitter` class that
  parses a PDB/mmCIF structure once, detects every disordered
  (ALTLOC-bearing) atom position and its altloc codes/occupancies, and
  splits the structure into one complete structure file per distinct
  altloc label found (`split()`), plus a JSON summary of what was found
  (`summary()`).
- Per-label atom resolution falls back to the highest-occupancy conformer
  (ties broken alphabetically) when a disordered position doesn't carry
  the requested label; every fallback is recorded in the JSON summary's
  `outputs.<label>.fallback_atoms` rather than silently dropping atoms.
- Each label's output is written from a fresh re-parse of the source file
  (not a shared in-memory structure) since Biopython's disordered-atom
  unpacking mutates the underlying atom objects during selection.
- Written altloc codes are blanked (`" "`/`.`) in every output — each
  file now represents one complete conformer, not a mixture; occupancy
  values are left untouched.
- `--format {pdb,cif}` output format, default `cif`.
- `histo-altloc-cleaner` CLI (Click-based, Rich console table output)
  with `--output-dir` and `--format` options.
- Claude Code / Claude Desktop skill (`skills/histo-altloc-cleaner/`)
  wrapping the CLI.
- Test suite (pytest) against a committed real structure with genuine
  ALTLOC A/B records (`7s79_1_aligned.cif`).
- `README.md`, `CLAUDE.md`, and design plan (`docs/PLAN.md`).
