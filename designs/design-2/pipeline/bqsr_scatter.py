"""Native chromosome-disjoint BQSR application and explicit record accounting.

The probe compares native whole-file and disjoint-interval application on actual
normal records, including placed and unplaced unmapped reads. It is an I/O and
quality-transformation check, not a new patient variant analysis.
"""
from pathlib import Path
import collections
import hashlib
import itertools
import json
import sys

import pysam
from common import gatk, run, receipt, save, sha


def record_accounting(path):
    with pysam.AlignmentFile(path, "rb") as bam:
        bam.check_index()
        return dict(contigs={x.contig: [x.mapped, x.unmapped]
                             for x in bam.get_index_statistics()},
                    nocoordinate=bam.nocoordinate,
                    samples=sorted({x["SM"] for x in bam.header.to_dict()["RG"]}))


def main():
    root=Path(sys.argv[1]);action=sys.argv[2]
    base=json.loads((root/"bqsr-plan.json").read_text())
    normal=Path(base["normal_bam"]);ref=Path(base["reference"])
    table=root/"bqsr/gather/normal.recalibration.table"
    done=json.loads((table.parent/"complete.json").read_text())
    assert done["status"]=="COMPLETED" and sha(table)==done["outputs"][table.name]["sha256"]
    if action=="probe":
        out=root/"bqsr/scatter-probe";out.mkdir(parents=True,exist_ok=True)
        scratch=Path("/tmp/bqsr-scatter-probe");scratch.mkdir(exist_ok=True)
        unsorted=scratch/"sample.unsorted.bam";sample=scratch/"sample.bam"
        categories=collections.Counter();seen=set()
        with pysam.AlignmentFile(normal,"rb") as source, pysam.AlignmentFile(unsorted,"wb",header=source.header) as dest:
            def add(read):
                key=read.to_string()
                if key in seen:return
                seen.add(key);dest.write(read)
                categories["unplaced_unmapped" if read.reference_id<0 else
                           "placed_unmapped" if read.is_unmapped else "mapped"]+=1
                if read.is_supplementary:categories["supplementary"]+=1
                if read.is_duplicate:categories["duplicate"]+=1
            # Every contig remains in the header; actual nonempty contigs are sampled.
            for chrom in source.references:
                for read in itertools.islice(source.fetch(chrom),100):add(read)
            source.reset()
            for read in itertools.islice(source.fetch(until_eof=True),1500000):
                if read.is_unmapped and read.reference_id>=0:add(read)
                if categories["placed_unmapped"]>=10:break
            for read in itertools.islice(source.fetch("*"),1000):add(read)
        assert categories["mapped"] and categories["placed_unmapped"] and categories["unplaced_unmapped"],categories
        run(out,"sample-sort",["samtools","sort","-o",sample,unsorted])
        run(out,"sample-index",["samtools","index",sample])
        whole=scratch/"whole.bam"
        run(out,"whole",gatk("ApplyBQSR","-R",ref,"-I",sample,"--bqsr-recal-file",table,"-O",whole,heap="3g"))
        with pysam.AlignmentFile(sample,"rb") as bam:
            chroms=list(bam.references)
        # Three disjoint dictionary-ordered groups plus the no-coordinate tail.
        groups=[chroms[:12],chroms[12:25],chroms[25:],["unmapped"]]
        parts=[]
        for index,group in enumerate(groups):
            if not group:continue
            part=scratch/f"part-{index}.bam";parts.append(part)
            run(out,f"part-{index}",gatk("ApplyBQSR","-R",ref,"-I",sample,
                 "--bqsr-recal-file",table,*[a for chrom in group for a in ["-L",chrom]],"-O",part,heap="3g"))
        combined=scratch/"gathered.bam"
        run(out,"gather",gatk("GatherBamFiles",*[a for p in parts for a in ["-I",p]],"-O",combined))
        def records(path):
            with pysam.AlignmentFile(path,"rb") as bam:
                return sorted(read.to_string() for read in bam.fetch(until_eof=True))
        raw=records(sample);a=records(whole);b=records(combined)
        assert len(raw)==len(a)==len(b) and a==b,"Native interval application did not preserve identical recalibrated records"
        # Alignment identity/sequence must stay unchanged. Native ApplyBQSR may
        # change quality-related optional tags; audit these separately. Exact
        # whole-versus-interval equality above includes QUAL and every tag.
        def identity(rows):
            return sorted(tuple(line.split("\t")[:10]) for line in rows)
        assert identity(raw)==identity(a)
        raw_ordered=sorted(raw,key=lambda s:s.split("\t")[:10])
        new_ordered=sorted(a,key=lambda s:s.split("\t")[:10])
        tag_changes=collections.Counter()
        for before,after in zip(raw_ordered,new_ordered):
            def tags(line):return {t.split(":",1)[0]:t for t in line.split("\t")[11:]}
            bt,at=tags(before),tags(after)
            for key in bt.keys()|at.keys():
                if bt.get(key)!=at.get(key):tag_changes[key]+=1
        save(out/"equivalence.json",dict(status="PASSED",scope="Actual normal-record native application equivalence; no variant inference",
             sampled_record_categories=dict(categories),record_count=len(raw),identical_recalibrated_records=True,
             original_record_identity_preserved=True,reference_sha256=base["reference_sha256"],
             native_whole_file_optional_tag_changes=dict(tag_changes),
             recalibration_table_sha256=sha(table),groups=groups,
             calibrated_record_multiset_sha256=hashlib.sha256("\n".join(a).encode()).hexdigest()))
        receipt(out,[out/"equivalence.json"],worker_sha256=sha(__file__))
    elif action=="apply":
        plan=json.loads((root/"bqsr-scatter-plan.json").read_text())
        probe=json.loads((root/"bqsr/scatter-probe/equivalence.json").read_text())
        assert probe["status"]=="PASSED" and probe["identical_recalibrated_records"]
        shard=next(x for x in plan["shards"] if x["id"]==sys.argv[3])
        out=root/"bqsr/scatter-apply"/shard["id"];out.mkdir(parents=True,exist_ok=True)
        bam=out/"recal.bam"
        run(out,"apply-bqsr",gatk("ApplyBQSR","-R",ref,"-I",normal,"--bqsr-recal-file",table,
             *[a for chrom in shard["intervals"] for a in ["-L",chrom]],"-O",bam,heap="6g"))
        run(out,"quickcheck",["samtools","quickcheck","-v",bam])
        run(out,"index",["samtools","index","-@","2",bam])
        save(out/"read-accounting.json",record_accounting(bam))
        receipt(out,[bam,Path(str(bam)+".bai"),out/"read-accounting.json"],shard=shard,
                plan_sha256=sha(root/"bqsr-scatter-plan.json"),worker_sha256=sha(__file__))
    elif action=="gather":
        plan=json.loads((root/"bqsr-scatter-plan.json").read_text())
        out=root/"bqsr/apply";out.mkdir(parents=True,exist_ok=True)
        assert not (out/"complete.json").exists(),"A completed calibrated BAM already exists"
        parts=[]
        for shard in plan["shards"]:
            directory=root/"bqsr/scatter-apply"/shard["id"]
            rec=json.loads((directory/"complete.json").read_text());assert rec["status"]=="COMPLETED"
            bam=directory/"recal.bam"
            assert bam.stat().st_size==rec["outputs"][bam.name]["bytes"]
            parts.append(bam)
        partial=out/"normal.recal.scattered.incomplete.bam";final=out/"normal.recal.bam"
        assert not final.exists()
        run(out,"gather-bams",gatk("GatherBamFiles",*[a for p in parts for a in ["-I",p]],"-O",partial,heap="8g"))
        run(out,"quickcheck",["samtools","quickcheck","-v",partial])
        run(out,"index",["samtools","index","-@","8",partial])
        original=record_accounting(normal);current=record_accounting(partial)
        assert original==current,"Disjoint BQSR application changed per-contig or unmapped record counts"
        save(out/"read-accounting.json",dict(before=original,after=current,identical=True))
        partial.rename(final);Path(str(partial)+".bai").rename(Path(str(final)+".bai"))
        result=receipt(out,[final,Path(str(final)+".bai"),out/"read-accounting.json"],
              normal_bam=str(final),normal_sample="TN25-112103-N",reference=str(ref),BQSR_applied=True,
              all_original_records_preserved=True,normal_original_sha256=base["normal_bam_sha256"],
              plan_sha256=sha(root/"bqsr-plan.json"),scatter_plan_sha256=sha(root/"bqsr-scatter-plan.json"),
              application="Native ApplyBQSR on disjoint whole-contig groups plus no-coordinate tail; native GatherBamFiles",
              equivalence_sha256=sha(root/"bqsr/scatter-probe/equivalence.json"),worker_sha256=sha(__file__))
        save(root/"bqsr/complete.json",result)
    else:raise ValueError(action)


if __name__=="__main__":main()
