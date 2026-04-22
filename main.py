"""
Process EMI archives: extract mzML files and run MZMine on each.

For each mzML file in the archives:
1. Extract it to a dedicated directory named after the file (without extension)
2. Run MZMine with the emi.mzbatch batch file
"""

from environments_utils import is_macos, is_linux
import subprocess
import sys
from pathlib import Path

from sirius_ms2_vs_ms12.constants import (
    MZMINE_LINUX_EXECCUTABLE,
    MZMINE_MACOS_EXECCUTABLE,
)

ARCHIVES: list[tuple[str, Path]] = [
    ("archive/emi-negative-archive.tar.zst", Path("negative")),
    ("archive/emi-positive-archive.tar.zst", Path("positive")),
]
BATCH_FILE: str = "src/sirius_ms2_vs_ms12/mzmine/emi.mzbatch"


def get_mzmine_executable() -> Path:
    """Return the MZMine executable path for the current OS."""
    if is_macos():
        return Path(MZMINE_MACOS_EXECCUTABLE)
    elif is_linux():
        return Path(MZMINE_LINUX_EXECCUTABLE)
    else:
        raise RuntimeError(f"Unsupported operating system.")


def list_mzml_members(archive_path: Path) -> list[str]:
    """Return all mzML member paths inside a .tar.zst archive."""
    result = subprocess.run(
        [
            "tar",
            "--use-compress-program=unzstd",
            "-tf",
            str(archive_path),
        ],
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
    Extract a single mzML member from a .tar.zst archive into dest_dir.

    The file is written to dest_dir/<stem>/<filename>.mzML, where <stem>
    is the filename without the .mzML extension.

    Returns the path to the extracted mzML file.
    """
    stem = Path(member).stem  # filename without .mzML
    name = Path(member).name  # filename with .mzML
    target_dir = dest_dir / stem
    target_dir.mkdir(parents=True, exist_ok=True)
    target_file = target_dir / name

    print(f"  Extracting {name} → {target_file}")

    result = subprocess.run(
        [
            "tar",
            "--use-compress-program=unzstd",
            "-xOf",
            str(archive_path),
            member,
        ],
        capture_output=True,
        check=True,
    )

    target_file.write_bytes(result.stdout)
    return target_file


def run_mzmine(
    mzmine_exe: Path,
    mzml_file: Path,
    output_dir: Path,
    batch_file: Path,
) -> None:
    """Run MZMine on a single mzML file."""
    cmd: list[str] = [
        str(mzmine_exe),
        "-i",
        str(mzml_file),
        "-o",
        str(output_dir),
        "-b",
        str(batch_file),
    ]
    print(f"  Running MZMine: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def process_archive(
    archive_path: Path,
    output_root: Path,
    mzmine_exe: Path,
    batch_file: Path,
) -> None:
    """Extract all mzML files from an archive and run MZMine on each."""
    print(f"\nProcessing archive: {archive_path}")

    members = list_mzml_members(archive_path)
    if not members:
        print("  No mzML files found — skipping.")
        return

    print(f"  Found {len(members)} mzML file(s).")

    for member in members:
        mzml_file = extract_mzml(archive_path, member, output_root)
        run_mzmine(
            mzmine_exe=mzmine_exe,
            mzml_file=mzml_file,
            output_dir=mzml_file.parent,
            batch_file=batch_file,
        )
        mzml_file.unlink()
        mgf_file = Path(mzml_file.stem + "_sirius.mgf")



def main() -> None:
    repo_root = Path(__file__).parent.resolve()

    mzmine_exe = get_mzmine_executable()
    if not mzmine_exe.exists():
        print(f"ERROR: MZMine executable not found at {mzmine_exe}", file=sys.stderr)
        sys.exit(1)

    batch_file = repo_root / BATCH_FILE
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
