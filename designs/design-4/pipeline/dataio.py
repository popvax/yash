"""Local, deterministic artifact IO; no network or patient-data uploads."""
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path


def resolve_path(path):
    p = Path(path)
    if not p.exists() and p.suffix != ".gz":
        gz = Path(str(p) + ".gz")
        if gz.exists():
            return gz
    return p


def resolve_gencode_fasta(filename, default_dir=None):
    """Locate an external GENCODE translation archive if present locally."""
    name = Path(filename).name
    candidates = []
    env_dir = os.environ.get("YASH_GENCODE_DIR")
    if env_dir:
        candidates.append(Path(env_dir) / name)
    if default_dir is not None:
        candidates.append((Path(default_dir) / filename).resolve())
    here = Path(__file__).resolve().parents[1]
    candidates.append(here / "references" / name)
    candidates.append(here.parents[1] / "work/gencode-references" / name)
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return candidates[1] if len(candidates) > 1 else candidates[0]


def sha(path):
    p = Path(path)
    if not p.exists() and p.suffix != ".gz":
        gz = Path(str(p) + ".gz")
        if gz.exists():
            return hashlib.sha256(gzip.decompress(gz.read_bytes())).hexdigest()
    digest = hashlib.sha256()
    with p.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_text_auto(path):
    p = resolve_path(path)
    if p.suffix == ".gz":
        return gzip.decompress(p.read_bytes()).decode("utf-8")
    return p.read_text()


def write_text_auto(path, text):
    p = resolve_path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    payload = text.encode("utf-8")
    if p.suffix == ".gz":
        p.write_bytes(gzip.compress(payload, compresslevel=9, mtime=0))
    else:
        p.write_bytes(payload)


def iter_rows(path):
    path = resolve_path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as stream:
        yield from csv.DictReader(
            (line for line in stream if not line.startswith("#")), delimiter="\t"
        )


def rows(path):
    return list(iter_rows(path))


def cell(value):
    if isinstance(value, (dict, list)):
        return json.dumps(value, sort_keys=True, separators=(",", ":"))
    return "" if value is None else str(value)


def table(path, records, fields=None):
    import io
    path = resolve_path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(dict.fromkeys(key for record in records for key in record))
    if not fields:
        raise ValueError("A table requires explicit fields or nonempty records")
    buf = io.StringIO(newline="")
    writer = csv.DictWriter(buf, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows({key: cell(record.get(key)) for key in fields} for record in records)
    payload = buf.getvalue().encode("utf-8")
    if path.suffix == ".gz":
        # mtime=0 removes timestamp drift in Python gzip.compress; note RFC 1952 OS byte (9) and deflate
        # bitstreams can still differ across external gzip/zlib implementations if regenerated on other stacks.
        path.write_bytes(gzip.compress(payload, compresslevel=9, mtime=0))
    else:
        path.write_bytes(payload)


def save(path, value):
    path = resolve_path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n").encode("utf-8")
    if path.suffix == ".gz":
        path.write_bytes(gzip.compress(payload, compresslevel=9, mtime=0))
    else:
        path.write_bytes(payload)


def fasta(path, label, sequence):
    Path(path).write_text(">" + label + "\n" + "".join(
        sequence[start:start + 80] + "\n" for start in range(0, len(sequence), 80)
    ))


def load(path):
    p = resolve_path(path)
    if p.suffix == ".gz":
        return json.loads(gzip.decompress(p.read_bytes()).decode("utf-8"))
    return json.loads(p.read_text())
