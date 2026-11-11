"""Check sample identity on the indexed normal BAM while calibration trains.

The final identity/contamination stage repeats fingerprinting after BQSR. This
early check is identity evidence only, not a somatic-variant assessment.
"""
from pathlib import Path
import csv
import json
import sys
import shutil

from common import gatk,run,receipt,save,sha
from fingerprint_panel import compatible_panel


def main():
    root=Path(sys.argv[1]);out=root/'identity_precal';out.mkdir(parents=True,exist_ok=True)
    previous=out/'fingerprint-precal.log'
    if previous.exists() and not (out/'complete.json').exists():
        archive=out/'attempts'/sha(previous)[:16];archive.mkdir(parents=True,exist_ok=True)
        for file in list(out.iterdir()):
            if file.is_file():shutil.copyfile(file,archive/file.name)
    alignment=json.loads((root/'alignment/complete.json').read_text())
    assert alignment['status']=='COMPLETED'
    assert json.loads((root/'resources/complete.json').read_text())['status']=='COMPLETED'
    normal=Path(alignment['normal_bam']);tumour=Path('/data/inputs/TN25-112103/DNA_TN25-112103.bam')
    metrics=out/'fingerprints.metrics.tsv';panel,panel_audit=compatible_panel(root,tumour,normal)
    run(out,'fingerprint-precal',gatk('CrosscheckFingerprints','-I',tumour,'-I',normal,'-H',panel,'-O',metrics,
        '--CROSSCHECK_BY','SAMPLE','--EXPECT_ALL_GROUPS_TO_MATCH','true','--LOD_THRESHOLD','-5',
        '--CALCULATE_TUMOR_AWARE_RESULTS','true','--REQUIRE_INDEX_FILES','true','--NUM_THREADS','2',
        '--VALIDATION_STRINGENCY','LENIENT'))
    lines=[line for line in metrics.read_text().splitlines() if line.strip() and not line.startswith('#')]
    rows=list(csv.DictReader(lines,delimiter='\t'))
    cross=[r for r in rows if {r['LEFT_GROUP_VALUE'],r['RIGHT_GROUP_VALUE']}=={'TN25-112103','TN25-112103-N'}]
    assert cross and all(r['RESULT']=='EXPECTED_MATCH' and float(r['LOD_SCORE'])>=5 for r in cross),cross
    save(out/'fingerprint-verdict.json',dict(status='MATCH',native_results=cross,pre_BQSR=True,
        scope='Identity evidence on the quality-checked duplicate-marked BAM; not proof of biological normality or somatic status'))
    receipt(out,[metrics,out/'fingerprint-verdict.json'],identity='MATCH',pre_BQSR=True,
            normal_bam_sha256=alignment['outputs'][normal.name]['sha256'],fingerprint_panel_sha256=sha(panel),
            fingerprint_coordinate_audit_sha256=sha(panel_audit),worker_sha256=sha(__file__))


if __name__=='__main__':main()
