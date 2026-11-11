"""Package the exact recorded source/evidence for historical stage reruns."""
from pathlib import Path
import gzip,io,json,tarfile
from yash_vaccine.provenance import RecordedSource


def package(repo,output):
    repo=Path(repo).resolve();output=Path(output);output.mkdir(parents=True,exist_ok=True)
    recorded=RecordedSource(repo)
    manifest_bytes=(repo/'designs/design-2/execution/frozen-inputs.json').read_bytes()
    manifest=json.loads(manifest_bytes)
    target=output/'frozen-inputs.json';archive=output/'frozen-inputs.tar.gz'
    payload={name:recorded.read(name,meta['sha256'],meta['bytes']) for name,meta in manifest['files'].items()}
    if target.exists() or archive.exists():
        assert target.exists() and archive.exists(),'Inspect incomplete recorded-input package'
        assert target.read_bytes()==manifest_bytes,'Use a new run root for different inputs'
        with tarfile.open(archive) as source:
            assert {m.name for m in source.getmembers()}==set(payload)
            for member in source.getmembers():assert member.isfile() and source.extractfile(member).read()==payload[member.name]
        return [target,archive]
    with archive.open('wb') as raw,gzip.GzipFile(filename='',mode='wb',mtime=0,fileobj=raw) as zipped,tarfile.open(fileobj=zipped,mode='w') as source:
        for name,body in payload.items():
            member=tarfile.TarInfo(name);member.size=len(body);member.mode=0o644;member.mtime=0
            source.addfile(member,io.BytesIO(body))
    target.write_bytes(manifest_bytes)
    return [target,archive]
