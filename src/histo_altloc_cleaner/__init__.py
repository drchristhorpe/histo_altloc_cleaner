"""Split a structure file with ALTLOC records into one file per conformer."""

from histo_altloc_cleaner.core import AltlocSplitter, StructureError, load_structure, split_altlocs

__all__ = ["AltlocSplitter", "StructureError", "load_structure", "split_altlocs"]
