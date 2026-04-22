"""
Process EMI archives: extract mzML files, run MZMine, then run Sirius
(via metabolite-annotator) in two modes for each generated MGF.

Pipeline per mzML:
  1. Extract mzML from archive  →  output/<stem>/<stem>.mzML
  2. Run MZMine                 →  output/<stem>/<stem>_sirius.mgf
  3. Run Sirius (normal)        →  output/<stem>/sirius_normal/   + results.csv
  4. Run Sirius (ms2_only)      →  output/<stem>/sirius_ms2_only/ + results.csv
"""


import subprocess
import sys
from enum import StrEnum
from pathlib import Path
from environments_utils import is_linux, is_macos
from sirius_ms2_vs_ms12.constants import (
    MZMINE_LINUX_EXECCUTABLE,
    MZMINE_MACOS_EXECCUTABLE,
)

ARCHIVES: list[tuple[str, Path]] = [
    ("archive/emi-positive-archive.tar.zst", Path("positive")),
    ("archive/emi-negative-archive.tar.zst", Path("negative")),
]
BATCH_FILE: str = "src/sirius_ms2_vs_ms12/mzmine/emi.mzbatch"
OUTPUT_ROOT: Path = Path("output")

# Sub-directory names for the two Sirius run modes
SIRIUS_NORMAL_DIR: str = "sirius_normal"
SIRIUS_MS2_ONLY_DIR: str = "sirius_ms2_only"

# CSV result filename produced inside each Sirius output directory
SIRIUS_RESULTS_CSV: str = "results.csv"


# ---------------------------------------------------------------------------
# Helpers: ion mode
# ---------------------------------------------------------------------------


class IonMode(StrEnum):
    """Ionization mode inferred from the archive filename."""

    POS = "pos"
    NEG = "neg"


def ion_mode_from_archive(archive_path: Path) -> IonMode:
    """
    Infer the ionization mode from the archive filename.

    Expects 'positive' or 'negative' to appear in the filename.
    Raises ValueError if neither is found.
    """
    name = archive_path.name.lower()
    if "positive" in name:
        return IonMode.POS
    if "negative" in name:
        return IonMode.NEG
    raise ValueError(
        f"Cannot determine ion mode from archive name: {archive_path.name}. "
        "Expected 'positive' or 'negative' in the filename."
    )


# ---------------------------------------------------------------------------
# Helpers: executables
# ---------------------------------------------------------------------------


def get_mzmine_executable() -> Path:
    """Return the MZMine executable path for the current OS."""
    if is_macos():
        return Path(MZMINE_MACOS_EXECCUTABLE)
    elif is_linux():
        return Path(MZMINE_LINUX_EXECCUTABLE)
    else:
        raise RuntimeError(f"Unsupported operating system.")


# ---------------------------------------------------------------------------
# Step 1 – archive extraction
# ---------------------------------------------------------------------------


def list_mzml_members(archive_path: Path) -> list[str]:
    """Return all mzML member paths present in a .tar.zst archive."""
    result = subprocess.run(
        ["tar", "--use-compress-program=unzstd", "-tf", str(archive_path)],
        capture_output=True,
        text=True,
        check=True,
    )
    return [
        line.strip()
        for line in result.stdout.splitlines()
        if line.strip().endswith(".mzML")
    ]


def extract_mzml(archive_path: Path, member: str, dest_dir: Path) -> Path:
    """
    Stream a single mzML member out of a .tar.zst archive.

    The file lands at ``dest_dir/<stem>/<filename>.mzML``.
    Returns the path to the written file.
    """
    stem: str = Path(member).stem
    name: str = Path(member).name
    target_dir: Path = dest_dir / stem
    target_dir.mkdir(parents=True, exist_ok=True)
    target_file: Path = target_dir / name

    print(f"  Extracting {name} → {target_file}")

    result = subprocess.run(
        ["tar", "--use-compress-program=unzstd", "-xOf", str(archive_path), member],
        capture_output=True,
        check=True,
    )
    target_file.write_bytes(result.stdout)
    return target_file


# ---------------------------------------------------------------------------
# Step 2 – MZMine
# ---------------------------------------------------------------------------


def run_mzmine(
    mzmine_exe: Path,
    mzml_file: Path,
    output_dir: Path,
    batch_file: Path,
) -> None:
    """Run MZMine on *mzml_file*, writing results into *output_dir*."""
    cmd: list[str] = [
        str(mzmine_exe),
        "-i", str(mzml_file),
        "-o", str(output_dir),
        "-b", str(batch_file),
    ]
    print(f"  Running MZMine: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def expected_mgf(mzml_file: Path) -> Path:
    """
    Return the MGF path that MZMine is expected to produce for *mzml_file*.

    Convention: ``<parent>/<stem>_sirius.mgf``
    """
    return mzml_file.parent / f"{mzml_file.stem}_sirius.mgf"


# ---------------------------------------------------------------------------
# Step 3 – Sirius (via metabolite-annotator)
# ---------------------------------------------------------------------------


def run_sirius(
    mgf_file: Path,
    output_dir: Path,
    ion_mode: IonMode,
    *,
    ms2_only: bool,
) -> Path:
    """
    Invoke ``uvx metabolite-annotator … sirius`` on *mgf_file*.

    Parameters
    ----------
    mgf_file:
        Path to the input MGF produced by MZMine.
    output_dir:
        Directory where Sirius will write its project space.
    ion_mode:
        Ionization mode (``pos`` or ``neg``).
    ms2_only:
        When *True*, pass ``--ms2_only`` to the sirius subcommand.

    Returns
    -------
    Path
        The *output_dir* that was passed to the tool (caller can inspect it).
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    cmd: list[str] = [
        "uvx",
        "metabolite-annotator",
        "-i", str(mgf_file),
        "-o", str(output_dir),
        "--ion-mode", ion_mode.value,
        "sirius",
    ]
    if ms2_only:
        cmd.append("--ms2-only")

    mode_label = "ms2_only" if ms2_only else "normal"
    print(f"  Running Sirius ({mode_label}): {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    return output_dir



# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------


def process_mzml(
    mzml_file: Path,
    mzmine_exe: Path,
    batch_file: Path,
    ion_mode: IonMode,
) -> None:
    """
    Full pipeline for a single mzML file:
      MZMine → MGF → Sirius (normal) → Sirius (ms2_only)
    """
    sample_dir: Path = mzml_file.parent

    # --- MZMine ---
    run_mzmine(
        mzmine_exe=mzmine_exe,
        mzml_file=mzml_file,
        output_dir=sample_dir,
        batch_file=batch_file,
    )

    mgf_file: Path = expected_mgf(mzml_file)
    if not mgf_file.exists():
        print(
            f"  WARNING: Expected MGF not found at {mgf_file}; skipping Sirius.",
            file=sys.stderr,
        )
        return

    # --- Sirius: normal run ---
    normal_dir: Path = sample_dir / SIRIUS_NORMAL_DIR
    run_sirius(mgf_file, normal_dir, ion_mode, ms2_only=False)

    # --- Sirius: ms2_only run ---
    ms2_only_dir: Path = sample_dir / SIRIUS_MS2_ONLY_DIR
    run_sirius(mgf_file, ms2_only_dir, ion_mode, ms2_only=True)


def process_archive(
    archive_path: Path,
    output_root: Path,
    mzmine_exe: Path,
    batch_file: Path,
) -> None:
    """Extract all mzML files from *archive_path* and run the full pipeline."""
    ion_mode: IonMode = ion_mode_from_archive(archive_path)
    print(f"\nProcessing archive: {archive_path}  (ion mode: {ion_mode.value})")

    members: list[str] = list_mzml_members(archive_path)
    if not members:
        print("  No mzML files found — skipping.")
        return

    print(f"  Found {len(members)} mzML file(s).")

    for member in members:
        mzml_file: Path = extract_mzml(archive_path, member, output_root)
        process_mzml(
            mzml_file=mzml_file,
            mzmine_exe=mzmine_exe,
            batch_file=batch_file,
            ion_mode=ion_mode,
        )


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main() -> None:
    repo_root: Path = Path(__file__).parent.resolve()

    mzmine_exe: Path = get_mzmine_executable()
    if not mzmine_exe.exists():
        print(f"ERROR: MZMine executable not found at {mzmine_exe}", file=sys.stderr)
        sys.exit(1)

    batch_file: Path = repo_root / BATCH_FILE
    if not batch_file.exists():
        print(f"ERROR: Batch file not found at {batch_file}", file=sys.stderr)
        sys.exit(1)

    for archive_str, output_subdir in ARCHIVES:
        output_root = repo_root / output_subdir
        output_root.mkdir(parents=True, exist_ok=True)
        archive_path = repo_root / archive_str
        if not archive_path.exists():
            print(
                f"WARNING: Archive not found, skipping: {archive_path}",
                file=sys.stderr,
            )
            continue
        process_archive(
            archive_path=archive_path,
            output_root=output_root,
            mzmine_exe=mzmine_exe,
            batch_file=batch_file,
        )

    print("\nDone.")


if __name__ == "__main__":
    main()