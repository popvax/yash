"""Audit current mutation targets with direct normal reads and native genotypes.

Normal-reference support alone is not labelled confirmed somatic. Final paired
calling is a separate required stage.
"""
from pathlib import Path
import collections
import csv
import json
import sys

import pysam
from common import gatk,run,receipt,save,sha

ROOT=Path(sys.argv[1]);OUT=ROOT/"audit";OUT.mkdir(parents=True,exist_ok=True)


def parse_metrics(target,bam):
    chrom,pos,ref,alt=target["allele"]
    assert len(ref)==len(alt)==1,"This specific 17-target audit expects single-base substitutions"
    fragments=collections.defaultdict(set)
    reads=collections.Counter()
    base_qualities=collections.defaultdict(list)
    for column in bam.pileup(chrom,pos-1,pos,truncate=True,stepper="nofilter",min_base_quality=0,
                             max_depth=1000000,ignore_overlaps=False,compute_baq=False):
        if column.reference_pos!=pos-1:continue
        for entry in column.pileups:
            read=entry.alignment
            if (read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate
                    or read.is_qcfail or read.mapping_quality<20):continue
            offset=entry.query_position
            if offset is None or read.query_qualities is None:continue
            quality=read.query_qualities[offset]
            if quality<20:continue
            base=read.query_sequence[offset].upper()
            label="ALT" if base==alt else "REF" if base==ref else "OTHER"
            fragments[read.query_name].add(label);reads[label]+=1;base_qualities[label].append(int(quality))
    counts=collections.Counter(next(iter(labels)) if len(labels)==1 else "DISAGREEING_MATES" for labels in fragments.values())
    return dict(qualifying_fragment_count=len(fragments),usable_fragment_count=sum(n for label,n in counts.items() if label!='DISAGREEING_MATES'),fragment_counts=dict(counts),read_counts=dict(reads),
                base_quality_medians={label:sorted(values)[len(values)//2] for label,values in base_qualities.items()})


def main():
    alignment=json.loads((ROOT/"alignment/complete.json").read_text())
    identity=json.loads((ROOT/"identity/complete.json").read_text())
    calibrated=json.loads((ROOT/"bqsr/complete.json").read_text())
    assert alignment["status"]==identity["status"]=="COMPLETED" and identity["identity"]=="MATCH"
    assert calibrated["status"]=="COMPLETED" and calibrated["BQSR_applied"]
    normal=Path(calibrated["normal_bam"]);reference=Path(alignment["reference"])
    targets=json.loads((ROOT/"inputs/current-mutation-targets.json").read_text())["targets"]
    assert len(targets)==17
    fasta=pysam.FastaFile(reference);order={name:i for i,name in enumerate(fasta.references)}
    header=pysam.VariantHeader();header.add_meta("fileformat","VCFv4.2")
    for name,length in zip(fasta.references,fasta.lengths):header.contigs.add(name,length=length)
    sites=OUT/"current-targets.vcf.gz"
    with pysam.VariantFile(sites,"wz",header=header) as writer:
        for target in sorted(targets,key=lambda t:(order[t["allele"][0]],t["allele"][1])):
            c,p,ref,alt=target["allele"]
            assert fasta.fetch(c,p-1,p-1+len(ref)).upper()==ref
            writer.write(writer.new_record(contig=c,start=p-1,alleles=(ref,alt),id=target["event_id"]))
    pysam.tabix_index(str(sites),preset="vcf",force=True)
    regions=ROOT/"inputs/current-mutation-targets.bed"
    genotypes=OUT/"normal-target-genotypes.vcf.gz"
    run(OUT,"genotype-current-targets",gatk("HaplotypeCaller","-R",reference,"-I",normal,
        "-L",regions,"--alleles",sites,"--output-mode","EMIT_ALL_CONFIDENT_SITES",
        "--native-pair-hmm-threads","2","-O",genotypes))
    calls={}
    with pysam.VariantFile(genotypes) as variants:
        assert list(variants.header.samples)==["TN25-112103-N"]
        for record in variants:
            calls[record.contig,record.pos]=dict(ref=record.ref,alts=list(record.alts or []),
                                               genotype=dict(record.samples["TN25-112103-N"].items()),
                                               native_FILTER=list(record.filter),native_QUAL=record.qual)
    audit=[]
    with pysam.AlignmentFile(normal,"rb") as bam:
        for target in targets:
            chrom,pos,ref,alt=target["allele"]
            measured=parse_metrics(target,bam)
            call=calls.get((chrom,pos));gt=call["genotype"].get("GT") if call else None
            gq=call["genotype"].get("GQ") if call else None
            indexes=[i+1 for i,a in enumerate(call["alts"]) if a==alt] if call and call["ref"]==ref else []
            normal_target_present=bool(gt and any(i in indexes for i in gt if i is not None))
            complete_gt=bool(gt and all(i is not None for i in gt))
            valid_substitution=bool(call and call['ref']==ref and complete_gt
                and all(len(call['ref'] if i==0 else call['alts'][i-1])==1 for i in gt))
            adequate=measured["usable_fragment_count"]>=20 and gq is not None and gq>=20 and complete_gt \
                and not (set(call['native_FILTER'])-{'PASS','.'}) and valid_substitution
            normal_alt=measured['fragment_counts'].get('ALT',0)
            if adequate and normal_target_present and normal_alt>0:
                verdict="TARGET_ALLELE_PRESENT_IN_NORMAL"
            elif adequate and not normal_target_present and normal_alt<3:
                verdict="NORMAL_REFERENCE_SUPPORTED_PENDING_PAIRED_CALL"
            else:
                verdict="UNRESOLVED_NORMAL_EVIDENCE"
            audit.append(dict(**target,**measured,normal_genotype=call,normal_target_allele_in_genotype=normal_target_present,
                              normal_callability_pass=adequate,normal_verdict=verdict,
                              normal_ALT_fragment_count=measured["fragment_counts"].get("ALT",0),
                              normal_REF_fragment_count=measured["fragment_counts"].get("REF",0),
                              somatic_classification="PENDING_FULL_PAIRED_DISCOVERY"))
    save(OUT/"current-mutation-audit.json",audit)
    fields=["antigen_id","gene","protein_change","allele","qualifying_fragment_count",
            "normal_ALT_fragment_count","normal_REF_fragment_count","normal_genotype","normal_verdict",
            "somatic_classification"]
    with (OUT/"current-mutation-audit.tsv").open("w") as stream:
        writer=csv.DictWriter(stream,fieldnames=fields,delimiter="\t",extrasaction="ignore");writer.writeheader()
        for row in audit:writer.writerow({k:json.dumps(v) if isinstance(v,(dict,list,tuple)) else v for k,v in row.items()})
    receipt(OUT,[sites,genotypes,OUT/"current-mutation-audit.json",OUT/"current-mutation-audit.tsv"],
        counts=dict(collections.Counter(row["normal_verdict"] for row in audit)),
        normal_bam=str(normal),normal_bam_sha256=calibrated["outputs"][normal.name]["sha256"],
        scope="17 current substitutions only; not nomination for full paired discovery",
        normal_callability_rule="At least 20 unambiguous qualifying fragment names and complete compatible native normal genotype GQ>=20; reference support also requires fewer than three normal ALT names",
        read_policy=dict(mapping_quality_min=20,base_quality_min=20,secondary=False,supplementary=False,
                         duplicates=False,qcfail=False,discordant_overlapping_mates="separate category"),
        pre_BQSR=False,identity_confirmed=True,worker_sha256=sha(__file__))


if __name__=="__main__":main()
