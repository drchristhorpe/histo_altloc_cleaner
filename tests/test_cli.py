from pathlib import Path

from click.testing import CliRunner

from histo_altloc_cleaner.cli import main

FIXTURE = Path(__file__).parent / "fixtures" / "7s79_1_aligned.cif"


def test_cli_writes_expected_files(tmp_path):
    runner = CliRunner()
    result = runner.invoke(main, [str(FIXTURE), "--output-dir", str(tmp_path)])
    assert result.exit_code == 0, result.output

    stem = FIXTURE.stem
    assert (tmp_path / f"{stem}_altlocA.cif").exists()
    assert (tmp_path / f"{stem}_altlocB.cif").exists()
    assert (tmp_path / f"{stem}_altlocs.json").exists()
    assert "2" in result.output


def test_cli_format_option(tmp_path):
    runner = CliRunner()
    result = runner.invoke(main, [str(FIXTURE), "--output-dir", str(tmp_path), "--format", "pdb"])
    assert result.exit_code == 0, result.output

    stem = FIXTURE.stem
    assert (tmp_path / f"{stem}_altlocA.pdb").exists()
    assert (tmp_path / f"{stem}_altlocB.pdb").exists()


def test_cli_missing_file_errors():
    runner = CliRunner()
    result = runner.invoke(main, ["does_not_exist.cif", "--output-dir", "out"])
    assert result.exit_code != 0
