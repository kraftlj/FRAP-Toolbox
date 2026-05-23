from __future__ import annotations

import argparse
import hashlib
import shutil
from pathlib import Path
from urllib.request import urlopen
from zipfile import ZipFile

ARCHIVE_URL = (
    "https://zenodo.org/records/20344310/files/"
    "frap-toolbox-legacy-user-guide-test-data.zip?download=1"
)
ARCHIVE_SHA256 = "3c17d0453bd325d3032859af75e90146008a7612a7bcc6fa9b92d90a9d046ddc"
ARCHIVE_NAME = "frap-toolbox-legacy-user-guide-test-data.zip"
CHUNK_SIZE = 1024 * 1024


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(CHUNK_SIZE), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _download(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with urlopen(url) as response, destination.open("wb") as handle:
        shutil.copyfileobj(response, handle, length=CHUNK_SIZE)


def _safe_extract(archive: Path, destination: Path) -> None:
    root = destination.resolve()
    with ZipFile(archive) as zf:
        for member in zf.infolist():
            target = (destination / member.filename).resolve()
            if root != target and root not in target.parents:
                raise ValueError(f"Archive member escapes destination: {member.filename}")
        zf.extractall(destination)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download, verify, and unpack the optional FRAP-Toolbox test-data archive."
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="Repository root where test-data/ should be restored.",
    )
    parser.add_argument(
        "--archive",
        type=Path,
        help="Use an existing archive instead of downloading from Zenodo.",
    )
    parser.add_argument(
        "--download-url",
        default=ARCHIVE_URL,
        help="Archive URL. Defaults to the published Zenodo record file.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Replace an existing test-data/ directory.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    repo_root = args.repo_root.resolve()
    test_data_dir = repo_root / "test-data"
    if test_data_dir.exists():
        if not args.force:
            raise SystemExit("test-data/ already exists. Use --force to replace it.")
        shutil.rmtree(test_data_dir)

    archive = args.archive
    if archive is None:
        archive = repo_root / "tmp" / ARCHIVE_NAME
        if not archive.exists():
            print(f"Downloading {args.download_url}")
            _download(args.download_url, archive)
    archive = archive.resolve()

    digest = _sha256(archive)
    if digest != ARCHIVE_SHA256:
        raise SystemExit(
            f"Checksum mismatch for {archive}: expected {ARCHIVE_SHA256}, got {digest}."
        )

    _safe_extract(archive, repo_root)
    if not test_data_dir.exists():
        raise SystemExit(f"Archive did not create {test_data_dir}.")
    print(f"Restored {test_data_dir}")


if __name__ == "__main__":
    main()
