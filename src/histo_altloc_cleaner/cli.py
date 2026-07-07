"""Command line interface for histo_altloc_cleaner."""

from __future__ import annotations

import click
from rich.console import Console
from rich.table import Table

from histo_altloc_cleaner.core import AltlocSplitter, StructureError


@click.command()
@click.argument("filename", type=click.Path(exists=True, dir_okay=False))
@click.option(
    "--output-dir",
    "-o",
    required=True,
    type=click.Path(file_okay=False, writable=True),
    help="Directory to write the split structure files and JSON summary into (created if missing).",
)
@click.option(
    "--format",
    "-f",
    "fmt",
    type=click.Choice(["pdb", "cif"]),
    default="cif",
    show_default=True,
    help="Output structure file format.",
)
def main(filename: str, output_dir: str, fmt: str) -> None:
    """Split a structure file with ALTLOC records into one complete
    structure file per alternate conformation.

    FILENAME is the path to a .cif/.mmcif or .pdb/.ent structure file
    containing ALTLOC records.
    """
    console = Console()

    try:
        splitter = AltlocSplitter(filename)
        result = splitter.split(output_dir, fmt=fmt)
    except StructureError as exc:
        raise click.ClickException(str(exc)) from exc

    summary = result["summary"]
    # Count disordered residues per label (a residue counts for a label if
    # any of its disordered atoms carries that label).
    residues_by_label: dict[str, int] = {label: 0 for label in summary["altloc_labels"]}
    for residue in summary["residues"]:
        labels_here = {label for atom in residue["disordered_atoms"].values() for label in atom}
        for label in labels_here:
            residues_by_label[label] += 1

    table = Table(title=f"ALTLOC split: {filename}")
    table.add_column("label")
    table.add_column("output file")
    table.add_column("disordered residues")
    table.add_column("fallback atoms")

    for label, path in result["labels"].items():
        fallback_count = len(summary["outputs"][label]["fallback_atoms"])
        table.add_row(label, str(path), str(residues_by_label[label]), str(fallback_count))

    console.print(table)
    console.print(
        f"[bold]{len(summary['altloc_labels'])}[/bold] altloc(s) "
        f"({', '.join(summary['altloc_labels'])}) written; summary at {result['summary_path']}"
    )


if __name__ == "__main__":
    main()
