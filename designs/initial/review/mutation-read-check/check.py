"""Supplementary read-only SNV evidence; this does not change GATK/Isovar gates."""
import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from statistics import median

import pysam


def eligible(read):
    return not (read.is_unmapped or read.is_secondary or read.is_supplementary
                or read.is_duplicate or read.is_qcfail) and read.mapping_quality >= 20 \
        and (not read.has_tag('NH') or read.get_tag('NH') == 1)


def inspect(bam, allele):
    chrom, position, ref, alt = allele
    assert len(ref) == len(alt) == 1, 'This supplementary check only handles SNVs'
    reads, fragments, qualities, orientations = Counter(), defaultdict(set), defaultdict(list), defaultdict(Counter)
    for read in bam.fetch(chrom, position - 1, position):
        if not eligible(read):
            continue
        for query_position, reference_position in read.get_aligned_pairs(matches_only=True):
            if reference_position != position - 1:
                continue
            quality = read.query_qualities[query_position] if read.query_qualities is not None else 0
            if quality < 20:
                continue
            base = read.query_sequence[query_position]
            reads[base] += 1
            fragments[base].add(read.query_name)
            qualities[base].append(quality)
            orientations[base]['reverse' if read.is_reverse else 'forward'] += 1
    return dict(read_counts=dict(reads), fragment_names_per_base={k: len(v) for k, v in fragments.items()},
                alt_BQ_median=median(qualities[alt]) if qualities[alt] else None,
                alt_strand=dict(orientations[alt]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--request', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    request = json.loads(args.request.read_text())
    result = dict(reviewed_utc=datetime.now(timezone.utc).isoformat(),
                  request_sha256=hashlib.sha256(args.request.read_bytes()).hexdigest(),
                  worker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  pysam_version=pysam.__version__,
                  policy='Primary nonduplicate QC-pass alignments; MQ >=20; NH=1 when present; base quality >=20. Supplementary observations, not new eligibility gates.',
                  mutation={})
    with pysam.AlignmentFile(request['DNA_BAM'], 'rb') as dna, pysam.AlignmentFile(request['RNA_BAM'], 'rb') as rna:
        for item in request['mutations']:
            result['mutation'][item['variant_id']] = dict(
                gene=item['gene'], allele=item['allele'], region_id=item['region_id'],
                DNA=inspect(dna, item['allele']), RNA=inspect(rna, item['allele']))
            print('CHECKED', item['gene'], flush=True)
    result['status'] = 'COMPLETED'
    args.output.write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
