"""Audit unknown platform metadata without copying or modifying the normal BAM."""
from pathlib import Path
import json
import sys

import pysam
from common import receipt,save,sha


def main():
    root=Path(sys.argv[1]);out=root/'normal_metadata';out.mkdir(parents=True,exist_ok=True)
    alignment=json.loads((root/'alignment/complete.json').read_text());assert alignment['status']=='COMPLETED'
    proof=root/'inputs/platform-default-verification.json';verification=json.loads(proof.read_text())
    assert verification['GATK_version']=='4.6.2.0'
    assert verification['identical_non_argument_recalibration_tables'] and verification['identical_recalibrated_read_records']
    assert verification['all_1000_records_preserved'] and not verification['normal_patient_calibration']
    assert all(r['training_exit']==0 and r['application_exit']==0 for r in verification['results'])
    assert all(r['returncode']==0 for r in verification['results'][1]['downstream_compatibility'].values())
    source=Path(alignment['normal_bam']);record=alignment['outputs'][source.name]
    assert source.stat().st_size==record['bytes']
    with pysam.AlignmentFile(source,'rb') as bam:
        bam.check_index();header=bam.header.to_dict()
        assert {group['SM'] for group in header['RG']}=={'TN25-112103-N'}
        assert all(not group.get('PL') for group in header['RG']),'Unexpected known platform; review rather than overriding it'
        counts=dict(contigs={x.contig:[x.mapped,x.unmapped] for x in bam.get_index_statistics()},nocoordinate=bam.nocoordinate)
    save(out/'metadata-audit.json',dict(original_header=header,original_index_counts=counts,header_modified=False,
         instrument_confirmation='User does not know the platform/model',recalibration_default_platform='UNKNOWN',
         native_argument=['--default-platform','UNKNOWN'],input_BAM_modified=False,
         source_alignment_SHA256=record['sha256'],verification_sha256=sha(proof),
         verification_scope=verification['scope'],new_patient_calibration_table_from_probe=False))
    receipt(out,[out/'metadata-audit.json'],normal_bam=str(source),normal_bam_bytes=record['bytes'],
        normal_bam_sha256=record['sha256'],source_alignment_sha256=record['sha256'],
        reference=alignment['reference'],reference_sha256=alignment['reference_sha256'],
        actual_instrument_unconfirmed=True,platform_metadata='UNKNOWN',original_header_unchanged=True,
        strategy='Native BaseRecalibrator --default-platform UNKNOWN; no separate header-copy BAM',worker_sha256=sha(__file__))


if __name__=='__main__':main()
