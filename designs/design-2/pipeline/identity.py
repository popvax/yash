"""Native tumour/normal fingerprinting, normal coverage and contamination QC."""
from pathlib import Path
import concurrent.futures
import csv
import json
import math
import sys

import pysam
from common import gatk,run,receipt,save,sha
from fingerprint_panel import compatible_panel

ROOT=Path(sys.argv[1]);OUT=ROOT/"identity";OUT.mkdir(parents=True,exist_ok=True)


def main():
    alignment=json.loads((ROOT/"alignment/complete.json").read_text())
    resources=json.loads((ROOT/"resources/complete.json").read_text())
    calibrated=json.loads((ROOT/"bqsr/complete.json").read_text())
    assert alignment["status"]==resources["status"]=="COMPLETED"
    assert calibrated["status"]=="COMPLETED" and calibrated["BQSR_applied"]
    normal=Path(calibrated["normal_bam"])
    tumour=Path("/data/inputs/TN25-112103/DNA_TN25-112103.bam")
    ref=Path(alignment["reference"])
    panel,panel_audit=compatible_panel(ROOT,tumour,normal)
    common=Path("/data/raft/references/homo_sapiens/vcfs/small_exac_common_3.hg38.vcf.gz")
    metrics=OUT/"fingerprints.metrics.tsv"
    run(OUT,"fingerprint",gatk("CrosscheckFingerprints","-I",tumour,"-I",normal,
        "-H",panel,"-O",metrics,"--CROSSCHECK_BY","SAMPLE",
        "--EXPECT_ALL_GROUPS_TO_MATCH","true","--LOD_THRESHOLD","-5",
        "--CALCULATE_TUMOR_AWARE_RESULTS","true","--REQUIRE_INDEX_FILES","true",
        "--NUM_THREADS","2","--VALIDATION_STRINGENCY","LENIENT"))
    lines=[line for line in metrics.read_text().splitlines() if line.strip() and not line.startswith("#")]
    rows=list(csv.DictReader(lines,delimiter="\t"))
    cross=[r for r in rows if {r["LEFT_GROUP_VALUE"],r["RIGHT_GROUP_VALUE"]}=={"TN25-112103","TN25-112103-N"}]
    assert cross,"Fingerprint metrics do not contain the intended tumour/normal comparison"
    assert all(r["RESULT"]=="EXPECTED_MATCH" for r in cross),"Sample identity needs review; do not start paired inference"
    assert all(float(r["LOD_SCORE"])>=5 for r in cross)
    early=ROOT/'identity_precal/complete.json'
    if early.exists():assert json.loads(early.read_text())['identity']=='MATCH','Resolve the preliminary identity result before proceeding'
    save(OUT/"fingerprint-verdict.json",dict(status="MATCH",native_results=cross,
         scope="Native haplotype identity evidence, not proof that either specimen is biologically normal"))
    pileups=OUT/"normal-pileups.table"
    run(OUT,"normal-pileups",gatk("GetPileupSummaries","-R",ref,"-I",normal,
        "-V",common,"-L",common,"-O",pileups))
    run(OUT,"normal-contamination",gatk("CalculateContamination","-I",pileups,
        "-O",OUT/"normal-contamination.table"))
    tumour_pileups=Path("/data/results/TN25-112103/mutation-evidence-v2-20261102/support/tumour-pileups.table")
    run(OUT,"matched-contamination",gatk("CalculateContamination","-I",tumour_pileups,
        "-matched",pileups,"-O",OUT/"tumour-contamination.table",
        "--tumor-segmentation",OUT/"allele-fraction-segments.table"))
    receipt(OUT,[metrics,OUT/"fingerprint-verdict.json",pileups,OUT/"normal-contamination.table",
                 OUT/"tumour-contamination.table",OUT/"allele-fraction-segments.table"],
        identity="MATCH",normal_bam=str(normal),normal_bam_sha256=calibrated["outputs"][normal.name]["sha256"],
        tumour_bam=str(tumour),normal_is_white_blood_cell_DNA="User supplied provenance",
        pre_BQSR=False,fingerprint_panel_sha256=sha(panel),fingerprint_coordinate_audit_sha256=sha(panel_audit),worker_sha256=sha(__file__))


if __name__=="__main__":main()
