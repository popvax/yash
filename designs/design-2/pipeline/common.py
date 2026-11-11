"""Small shared helpers; importing this module starts no work."""
from pathlib import Path
import datetime as dt
import hashlib
import json
import subprocess


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b""):
            digest.update(block)
    return digest.hexdigest()


def save(path, value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_suffix(path.suffix+".tmp")
    temporary.write_text(json.dumps(value,indent=2,default=str)+"\n")
    temporary.replace(path)


def run(directory, name, command, *, stdout=None):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    command=list(map(str,command))
    save(directory/(name+".command.json"),dict(command=command))
    save(directory/"status.json",dict(phase=name,updated_utc=now()))
    print(json.dumps(dict(phase=name,updated_utc=now())),flush=True)
    with (directory/(name+".log")).open("w") as log:
        subprocess.run(command,stdout=stdout or log,stderr=log,check=True)


def receipt(directory, paths, **metadata):
    directory=Path(directory)
    result=dict(status="COMPLETED",completed_utc=now(),**metadata,
                outputs={str(p.relative_to(directory)):dict(bytes=p.stat().st_size,sha256=sha(p)) for p in paths})
    save(directory/"complete.json",result)
    save(directory/"status.json",result)
    print(json.dumps(result),flush=True)
    return result


def gatk(tool,*args,heap="6g"):
    return ["/gatk/gatk","--java-options",f"-Xmx{heap} -XX:ParallelGCThreads=2 -XX:ConcGCThreads=1",tool,*map(str,args)]
