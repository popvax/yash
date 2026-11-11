"""Native GATK BQSR training/gather/application; no work at import."""
from pathlib import Path
import json
import sys

import pysam
from common import gatk,run,receipt,save,sha


def main():
    root=Path(sys.argv[1]);action=sys.argv[2]
    plan=json.loads((root/"bqsr-plan.json").read_text())
    normal=Path(plan["normal_bam"]);reference=Path(plan["reference"])
    assert normal.stat().st_size==plan["normal_bam_bytes"]
    assert json.loads((root/"alignment/complete.json").read_text())["status"]=="COMPLETED"
    assert json.loads((root/"resources/complete.json").read_text())["status"]=="COMPLETED"
    if action=="train":
        shard=next(x for x in plan["shards"] if x["id"]==sys.argv[3])
        out=root/"bqsr/training"/shard["id"];out.mkdir(parents=True,exist_ok=True)
        table=out/"recalibration.table"
        command=gatk("BaseRecalibrator","-R",reference,"-I",normal,"-L",shard["interval"],
                     "--default-platform",plan["recalibration_default_platform"],
                     "--known-sites",plan["known_sites"][0],"--known-sites",plan["known_sites"][1],"-O",table)
        run(out,"base-recalibrator",command)
        receipt(out,[table],shard=shard,normal_bam_sha256=plan["normal_bam_sha256"],
                plan_sha256=sha(root/"bqsr-plan.json"),recalibration_default_platform=plan["recalibration_default_platform"],worker_sha256=sha(__file__))
    elif action=="gather":
        out=root/"bqsr/gather";out.mkdir(parents=True,exist_ok=True)
        tables=[]
        for shard in plan["shards"]:
            directory=root/"bqsr/training"/shard["id"]
            rec=json.loads((directory/"complete.json").read_text());assert rec["status"]=="COMPLETED"
            table=directory/"recalibration.table"
            assert sha(table)==rec["outputs"][table.name]["sha256"]
            tables.append(table)
        final=out/"normal.recalibration.table"
        run(out,"gather-bqsr",gatk("GatherBQSRReports",*[arg for path in tables for arg in ["-I",path]],"-O",final,heap="10g"))
        receipt(out,[final],shards=len(tables),plan_sha256=sha(root/"bqsr-plan.json"),worker_sha256=sha(__file__))
    elif action=="apply":
        out=root/"bqsr/apply";out.mkdir(parents=True,exist_ok=True)
        training=root/"bqsr/gather/normal.recalibration.table"
        assert json.loads((training.parent/"complete.json").read_text())["status"]=="COMPLETED"
        partial=out/"normal.recal.incomplete.bam";final=out/"normal.recal.bam"
        run(out,"apply-bqsr",gatk("ApplyBQSR","-R",reference,"-I",normal,"--bqsr-recal-file",training,"-O",partial,heap="12g"))
        run(out,"quickcheck",["samtools","quickcheck","-v",partial])
        run(out,"index",["samtools","index","-@","4",partial])
        def counts(path):
            with pysam.AlignmentFile(path,"rb") as bam:
                bam.check_index()
                return dict(contigs={x.contig:[x.mapped,x.unmapped] for x in bam.get_index_statistics()},
                            nocoordinate=bam.nocoordinate,samples=sorted({x["SM"] for x in bam.header.to_dict()["RG"]}))
        original_counts=counts(normal);new_counts=counts(partial)
        assert original_counts==new_counts,"BQSR changed per-contig alignment/unmapped record counts"
        partial.rename(final)
        Path(str(partial)+".bai").rename(Path(str(final)+".bai"))
        save(out/"read-accounting.json",dict(before=original_counts,after=new_counts,identical=True))
        rec=receipt(out,[final,Path(str(final)+".bai"),out/"read-accounting.json"],
                    normal_bam=str(final),normal_sample="TN25-112103-N",reference=str(reference),
                    BQSR_applied=True,all_original_records_preserved=True,normal_original_sha256=plan["normal_bam_sha256"],
                    plan_sha256=sha(root/"bqsr-plan.json"),worker_sha256=sha(__file__))
        save(root/"bqsr/complete.json",rec)
    else:
        raise ValueError(action)


if __name__=="__main__":main()
