"""Build a checksum-listed local review ZIP; no uploads, model weights or credentials."""
import argparse
import re
import tempfile
import zipfile
from pathlib import Path

from dataio import save, sha

HERE = Path(__file__).resolve().parents[1]


def included_files():
    return [
        path for path in sorted(HERE.rglob("*"))
        if path.is_file() and path.name not in (
            "manifest.json", ".DS_Store",
            "gencode.v46.pc_translations.fa.gz", "gencode.v50.pc_translations.fa.gz",
        )
        and not any(part in ("__pycache__", "consensus-search") for part in path.relative_to(HERE).parts)
    ]


def write_manifest():
    files = included_files()
    pattern = re.compile(rb"AIza[A-Za-z0-9_-]{20,}|gh[pousr]_[A-Za-z0-9]{30,}")
    for path in files:
        if path.is_symlink():
            raise ValueError("Do not package symlinked external inputs: " + str(path.relative_to(HERE)))
        with path.open("rb") as stream:
            tail = b""
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                data = tail + block
                if pattern.search(data):
                    raise ValueError("Possible credential found; review file without exposing it: " + str(path.relative_to(HERE)))
                tail = data[-256:]
    save(HERE / "manifest.json", {
        "design_id": "design-4", "revision": "2026-11-15-source-and-reporting-audit",
        "scope": "Research artifacts only, not a clinically validated or complete therapeutic molecule.",
        "clinical_release_ready": False,
        "files_sha256": {str(path.relative_to(HERE)): sha(path) for path in files},
        "patient_data_notice": "Contains patient-derived genomic evidence from the source project's public release and supplied research artifacts. Respect the source data permissions and third-party terms; this program performs no upload.",
    })
    return files


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--manifest-only", action="store_true")
    args = parser.parse_args()
    files = write_manifest()
    if args.manifest_only:
        print("Wrote credential-screened manifest for", len(files), "files")
        return
    if args.output is None:
        parser.error("--output is required unless --manifest-only is used")
    destination = args.output.resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        raise ValueError("Use a new ZIP path; preserve any existing archive explicitly first")
    with tempfile.NamedTemporaryFile(prefix=".design4-review-", suffix=".zip",
                                     dir=destination.parent, delete=False) as handle:
        temporary = Path(handle.name)
    try:
        with zipfile.ZipFile(temporary, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
            for path in files + [HERE / "manifest.json"]:
                archive.write(
                    path, "designs/design-4/" + str(path.relative_to(HERE)),
                    compress_type=zipfile.ZIP_STORED if path.suffix == ".gz" else zipfile.ZIP_DEFLATED,
                )
        with zipfile.ZipFile(temporary) as archive:
            if archive.testzip() is not None:
                raise ValueError("ZIP CRC validation failed")
        temporary.chmod(0o600)
        temporary.replace(destination)
    finally:
        if temporary.exists():
            temporary.unlink()
    print("Local review ZIP:", destination, "sha256:", sha(destination))


if __name__ == "__main__":
    main()
