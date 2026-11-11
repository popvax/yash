"""Raw normal FASTQ -> untrimmed QC -> BWA-MEM2 -> fixmate/sort/markdup -> QC.

Runs only as an explicit Modal worker; it never changes the supplied inputs.
"""
from pathlib import Path
import datetime as dt
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
import threading
import time

import pysam

BWA = "/opt/bwa-mem2-2.3_x64-linux/bwa-mem2"
ROOT = Path(sys.argv[1])
D = ROOT / "alignment"
D.mkdir(parents=True, exist_ok=True)
START = time.monotonic()


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def save(name, data):
    target = D / name
    tmp = target.with_suffix(target.suffix + ".tmp")
    tmp.write_text(json.dumps(data, indent=2) + "\n")
    tmp.replace(target)


def status(phase, **extra):
    data = dict(phase=phase, updated_utc=now(), elapsed_seconds=time.monotonic()-START, **extra)
    save("status.json", data)
    print(json.dumps(data), flush=True)


def run(name, args, stdout=None):
    save(name + ".command.json", dict(command=list(map(str, args))))
    status(name)
    with (D / (name + ".log")).open("w") as log:
        subprocess.run(list(map(str, args)), stdout=stdout or log, stderr=log, check=True)


def main():
    if (D / "complete.json").exists():
        raise RuntimeError("Completed alignment already exists; verify it without overwriting")
    source = Path("/data/inputs/TN25-112103/matched-normal-haystack-20261107")
    intake = json.loads((source / "transfer-receipt.json").read_text())
    assert intake["status"] == "complete" and len(intake["objects"]) == 2
    paths = [Path(obj["destination"]) for obj in sorted(intake["objects"], key=lambda x:x["key"])]
    for obj in intake["objects"]:
        assert obj["checksum_verified"] and Path(obj["destination"]).stat().st_size == obj["size"]
    ref = Path("/data/results/TN25-112103/mutect2-tumour-only-20261029/resources/reference.fa")
    expected = "c3e61ad6d770db22a6f27fdbcd3dba93d0d80f649018c8887ec803f12ea90450"
    status("reference-verification")
    assert sha(ref) == expected
    with gzip.open(paths[0], "rt") as stream:
        first_header = stream.readline().strip()
    # Seven colon-separated fields are the standard Illumina header layout.
    # Keep PL unspecified if the available header does not establish it.
    platform = "ILLUMINA" if len(first_header.split()[0].split(":")) == 7 else None
    rg = "@RG\\tID:BatchEXT867\\tSM:TN25-112103-N\\tLB:BatchEXT867\\tPU:BatchEXT867"
    if platform:
        rg += "\\tPL:" + platform
    versions = {}
    for name, command in {"bwa_mem2":[BWA,"version"], "samtools":["samtools","--version"], "fastp":["fastp","--version"]}.items():
        p = subprocess.run(command, capture_output=True, text=True, errors="replace")
        versions[name] = (p.stdout + p.stderr)[:3000]
    save("inputs.json", dict(normal_inputs=intake["objects"], reference=str(ref),
         reference_sha256=expected, normal_sample="TN25-112103-N", read_group=rg,
         platform_inferred_from_header=platform,
         source_reference_caveat="Original Caris custom masks are unavailable; compatible GRCh38 analysis reference retained",
         software=versions, worker_sha256=sha(__file__)))
    index = ROOT / "resources/bwa-mem2/reference"
    index.parent.mkdir(parents=True, exist_ok=True)
    if not (D / "index.complete.json").exists():
        run("reference-index", [BWA, "index", "-p", index, ref])
        index_files = [{"path":str(p),"bytes":p.stat().st_size,"sha256":sha(p)}
                       for p in sorted(index.parent.glob("reference*")) if p.is_file()]
        assert index_files
        save("index.complete.json", dict(reference_sha256=expected, files=index_files))
    tmp = Path("/tmp/yash-normal-alignment")
    tmp.mkdir(exist_ok=True)
    sorted_bam = D / "normal.sorted.bam"
    normal = D / "normal.dedup.bam"
    if not sorted_bam.exists() and not normal.exists():
        commands = [
            ["fastp","-i",str(paths[0]),"-I",str(paths[1]),"--stdout",
             "--disable_adapter_trimming","--disable_quality_filtering",
             "--disable_length_filtering","--disable_trim_poly_g",
             "--thread","4","--json",str(D/"fastp.json"),"--html",str(D/"fastp.html")],
            [BWA,"mem","-p","-t","24","-K","100000000","-R",rg,str(index),"-"],
            ["samtools","view","-u","-"],
            ["samtools","fixmate","-m","-u","-","-"],
            ["samtools","sort","-@","6","-m","2G","-T",str(tmp/"sort"),"-o",str(sorted_bam),"-"],
        ]
        save("alignment-pipeline.command.json", dict(commands=commands,
             read_modification="None: trimming, correction, merging and filtering are disabled",
             deterministic_bwa_batch_bases=100000000))
        status("alignment")
        handles=[];processes=[]
        try:
            previous=None
            for i, command in enumerate(commands):
                log=(D/f"pipeline-{i}.log").open("w");handles.append(log)
                process=subprocess.Popen(command,stdin=previous,stdout=subprocess.PIPE if i<len(commands)-1 else log,stderr=log)
                if previous is not None:previous.close()
                previous=process.stdout;processes.append(process)
            while any(process.poll() is None for process in processes):
                failures=[(i,p.returncode) for i,p in enumerate(processes) if p.returncode not in [None,0]]
                if failures:raise RuntimeError("Alignment pipeline failed: "+repr(failures))
                status("alignment", output_bam_bytes=sorted_bam.stat().st_size if sorted_bam.exists() else 0,
                       progress_log_tail=(D/"pipeline-1.log").read_text(errors="replace")[-1000:])
                time.sleep(30)
            codes=[p.wait() for p in processes]
            assert all(code==0 for code in codes),codes
        finally:
            for p in processes:
                if p.poll() is None:p.terminate()
            for handle in handles:handle.close()
        run("sorted-quickcheck",["samtools","quickcheck","-v",sorted_bam])
        fq=json.loads((D/"fastp.json").read_text())
        before=fq["summary"]["before_filtering"];after=fq["summary"]["after_filtering"]
        assert before["total_reads"] == after["total_reads"] and before["total_bases"] == after["total_bases"]
        save("fastq-qc-verified.json", dict(total_reads=before["total_reads"],total_bases=before["total_bases"],
             q30_rate=before["q30_rate"],gc_content=before["gc_content"],read_count_preserved=True,
             base_count_preserved=True))
    if not normal.exists():
        run("mark-duplicates",["samtools","markdup","-@","8","-s","-f",D/"duplicate-metrics.txt",sorted_bam,normal])
    run("normal-quickcheck",["samtools","quickcheck","-v",normal])
    run("normal-index",["samtools","index","-@","8",normal])
    for name, command in {
        "flagstat":["samtools","flagstat","-@","8",normal],
        "stats":["samtools","stats","-@","8",normal],
        "idxstats":["samtools","idxstats",normal],
    }.items():
        with (D/(name+".txt")).open("w") as target:run(name,command,stdout=target)
    with pysam.AlignmentFile(normal,"rb") as bam:
        bam.check_index()
        with pysam.FastaFile(ref) as genome:
            assert all(genome.get_reference_length(c)==n for c,n in zip(bam.references,bam.lengths))
        sample_ids=sorted({rg["SM"] for rg in bam.header.to_dict().get("RG",[])})
        assert sample_ids == ["TN25-112103-N"]
        mapped=sum(x.mapped for x in bam.get_index_statistics())
    status("hashing-normal-bam")
    outputs={p.name:dict(bytes=p.stat().st_size,sha256=sha(p)) for p in [normal,Path(str(normal)+".bai"),D/"fastp.json",D/"fastq-qc-verified.json"]}
    save("complete.json",dict(status="COMPLETED",completed_utc=now(),normal_bam=str(normal),
         normal_sample="TN25-112103-N",reference=str(ref),reference_sha256=expected,
         mapped_alignment_records=mapped,outputs=outputs,
         identity_confirmed=False,BQSR_applied=False,
         next="Identity/contamination, coverage and paired variant calling; assess BQSR separately"))
    status("COMPLETED")


if __name__=="__main__":
    try:main()
    except Exception as exc:
        status("FAILED",error_type=type(exc).__name__,error=str(exc))
        raise
