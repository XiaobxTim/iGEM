from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import shutil
import tarfile
import tempfile
from urllib.request import urlopen


PHYSICELL_VERSION = "1.14.2"
SOURCE_URL = "https://github.com/MathCancer/PhysiCell/archive/refs/tags/1.14.2.tar.gz"
SOURCE_SHA256 = "de7ece9a990e9dc1d035a57b40c1a205e1b4a43155654577b5acad819f61894a"


def _digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def setup_physicell(destination: str | Path, archive: str | Path | None = None) -> Path:
    """Install the checksum-pinned PhysiCell source without replacing existing data."""

    destination = Path(destination).resolve()
    root = destination / "PhysiCell"
    version_file = root / "VERSION.txt"
    if version_file.exists() and version_file.read_text(encoding="utf-8").strip() == PHYSICELL_VERSION:
        return root
    if destination.exists():
        raise FileExistsError(f"refusing to replace incomplete destination: {destination}")

    destination.parent.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="physicell-setup-", dir=destination.parent))
    try:
        bundle = work / "PhysiCell.tar.gz"
        if archive is None:
            with urlopen(SOURCE_URL) as response, bundle.open("wb") as output:
                shutil.copyfileobj(response, output)
        else:
            shutil.copy2(Path(archive), bundle)
        actual = _digest(bundle)
        if actual != SOURCE_SHA256:
            raise ValueError(f"PhysiCell archive checksum mismatch: {actual}")
        unpacked = work / "unpacked"
        unpacked.mkdir()
        with tarfile.open(bundle, "r:gz") as tar:
            tar.extractall(unpacked, filter="data")
        source = unpacked / f"PhysiCell-{PHYSICELL_VERSION}"
        if not (source / "VERSION.txt").is_file():
            raise ValueError("PhysiCell archive has an unexpected layout")
        destination.mkdir()
        shutil.move(str(source), str(root))
        return root
    except Exception:
        if destination.exists():
            shutil.rmtree(destination)
        raise
    finally:
        shutil.rmtree(work, ignore_errors=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Install pinned PhysiCell 1.14.2")
    parser.add_argument(
        "--destination", type=Path, default=Path("/private/tmp/PhysiCell-1.14.2")
    )
    parser.add_argument("--archive", type=Path)
    arguments = parser.parse_args()
    print(setup_physicell(arguments.destination, arguments.archive))


if __name__ == "__main__":
    main()
