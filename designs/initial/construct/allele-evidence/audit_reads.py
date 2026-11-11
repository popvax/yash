"""Read-only SEC/MITD evidence audit using previously extracted HLA fragments.

Run on the existing Modal volume. This does not change a vaccine sequence.
Coordinates in per-codon output are 1-based allele-reference coordinates.
"""
import argparse
import collections
import csv
import hashlib
import json
import subprocess
from pathlib import Path

import Bio
from Bio import SeqIO
from Bio.Seq import Seq
import pysam

SOURCE = Path('/data/results/TN25-112103/hla-class2-20261028')
OUT = Path('/data/results/TN25-112103/trafficking-domain-audit-20261102')
OUT.mkdir(parents=True, exist_ok=True)
MIN_BQ = 30
MIN_MQ = 20
TARGETS = {
    'SEC_C0102': dict(allele='HLA-C*01:02:01:23', first=1, last=26,
                     protein='MRVMAPRTLILLLSGALALTETWACS'),
    'MITD_B1508': dict(allele='HLA-B*15:08:01:01', first=308, last=362,
                      protein='IVGIVAGLAVLAVVVIGAVVATVMCRRKSSGGKGGSYSQAASSDSAQGSDVSLTA'),
}


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def family(allele):
    return ':'.join(allele.split(':')[:2])


def read_families():
    # DNA provides representatives; add RNA-only families as competing decoys.
    choices = {}
    for source, name in [('dna', 'TN25-112103-DNA'), ('rna', 'TN25-112103')]:
        for line in (SOURCE / source / f'{name}_allele.tsv').read_text().splitlines():
            allele, quality = line.split()
            choices.setdefault(family(allele), allele)
    for target in TARGETS.values():
        choices[family(target['allele'])] = target['allele']
    return choices


def prepare_reference(kind, choices):
    source = SOURCE / 'hla-reference' / f'hla_{kind}_seq.fa'
    selected = {}
    for record in SeqIO.parse(source, 'fasta'):
        if record.id in choices.values():
            selected[record.id] = record
    assert set(choices.values()) <= selected.keys(), set(choices.values()) - selected.keys()
    # Exact duplicate references are collapsed, with all aliases retained. Otherwise
    # identical references would make an informative fragment appear ambiguous.
    groups = collections.defaultdict(list)
    for allele, record in selected.items():
        groups[str(record.seq)].append(allele)
    aliases, records = {}, []
    for sequence, names in groups.items():
        preferred = next((t['allele'] for t in TARGETS.values() if t['allele'] in names), names[0])
        record = selected[preferred]
        aliases[preferred] = names
        records.append(record)
    path = OUT / f'{kind}-competitive-reference.fasta'
    SeqIO.write(records, path, 'fasta')
    details = {}
    for name, target in TARGETS.items():
        rec = selected[target['allele']]
        fields = rec.description.split()
        n = int(fields[1])
        bounds = list(map(int, fields[2:]))
        assert len(bounds) == 2 * n
        # T1K FASTA headers encode inclusive 0-based coding-exon bounds.
        cds_positions = [p for i in range(n) for p in range(bounds[2*i], bounds[2*i+1]+1)]
        coding = ''.join(str(rec.seq)[p] for p in cds_positions)
        translated = str(Seq(coding).translate())
        first, last = target['first'], target['last']
        assert translated[first-1:last] == target['protein'], (kind, name, translated)
        details[name] = dict(**target, coding_sequence=coding,
                             coding_positions_0based=cds_positions,
                             reference_sequence=str(rec.seq),
                             reference_length=len(rec))
    return path, details, aliases


def run(cmd, **kwargs):
    print('RUN', ' '.join(map(str, cmd)), flush=True)
    subprocess.run(cmd, check=True, **kwargs)


def primary_groups(path):
    """Yield full fragments from a name-sorted BAM; secondary placements excluded."""
    with pysam.AlignmentFile(path, 'rb') as bam:
        previous, records = None, []
        for read in bam:
            if read.query_name != previous and records:
                yield previous, records
                records = []
            previous = read.query_name
            if not (read.is_secondary or read.is_supplementary or read.is_qcfail):
                records.append(read)
        if records:
            yield previous, records


def merged_bases(reads):
    """Combine mates; conflicting high-quality bases are excluded, not voted."""
    values = collections.defaultdict(set)
    for read in reads:
        for query_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
            if read.query_qualities[query_pos] >= MIN_BQ:
                values[ref_pos].add(read.query_sequence[query_pos])
    return {p: next(iter(b)) for p, b in values.items() if len(b) == 1}


def analyse(kind, path, targets):
    counts = {name: {
        'codons': [collections.Counter() for _ in t['protein']],
        'strict_codons': [collections.Counter() for _ in t['protein']],
        'starts': [collections.defaultdict(set) for _ in t['protein']],
        'bases': [collections.Counter() for _ in range(3*len(t['protein']))],
        'qualifying_fragments': 0,
        'full_domain_exact_fragments': 0,
        'distinct_full_domain_fragment_coordinates': set(),
        'phase_329_349_exact_fragments': 0,
        'strict_read_329_349_exact': 0,
        'dna_coordinate_collapsed_codons': [collections.defaultdict(set) for _ in t['protein']],
    } for name, t in targets.items()}
    totals = collections.Counter()
    for name, reads in primary_groups(path):
        totals['candidate_fragments'] += 1
        mapped = [r for r in reads if not r.is_unmapped]
        if not mapped:
            continue
        refs = {r.reference_name for r in mapped}
        # Both mates, if mapped, must agree on the allele reference, face opposite
        # directions, and include at least one unambiguous alignment. BWA can assign
        # an ambiguous mate through its informative partner.
        if len(refs) != 1 or max(r.mapping_quality for r in mapped) < MIN_MQ:
            continue
        if len(mapped) > 1 and len({r.is_reverse for r in mapped}) != 2:
            continue
        ref = next(iter(refs))
        coordinates = tuple(sorted((r.reference_start, r.reference_end, r.is_reverse) for r in mapped))
        totals['confident_reference_fragments'] += 1
        for target_name, target in targets.items():
            if ref != target['allele']:
                continue
            c = counts[target_name]
            c['qualifying_fragments'] += 1
            bases = merged_bases(mapped)
            strict_bases = merged_bases([r for r in mapped if r.mapping_quality >= MIN_MQ])
            all_positions = target['coding_positions_0based']
            lo, hi = 3*(target['first']-1), 3*target['last']
            positions = all_positions[lo:hi]
            expected = target['coding_sequence'][lo:hi]
            for i, pos in enumerate(positions):
                if pos in bases:
                    c['bases'][i][bases[pos]] += 1
            for i in range(len(target['protein'])):
                positions3 = positions[3*i:3*i+3]
                if all(p in strict_bases for p in positions3):
                    codon = ''.join(strict_bases[p] for p in positions3)
                    c['strict_codons'][i][str(Seq(codon).translate())] += 1
                if all(p in bases for p in positions3):
                    codon = ''.join(bases[p] for p in positions3)
                    aa = str(Seq(codon).translate())
                    c['codons'][i][aa] += 1
                    c['starts'][i][aa].add(coordinates)
                    c['dna_coordinate_collapsed_codons'][i][aa].add(coordinates)
            if all(bases.get(pos) == base for pos, base in zip(positions, expected)):
                c['full_domain_exact_fragments'] += 1
                c['distinct_full_domain_fragment_coordinates'].add(coordinates)
            if target_name == 'MITD_B1508':
                phase_positions = all_positions[3*328:3*349]
                phase_expected = target['coding_sequence'][3*328:3*349]
                if all(bases.get(pos) == base for pos, base in zip(phase_positions, phase_expected)):
                    c['phase_329_349_exact_fragments'] += 1
                for read in mapped:
                    if read.mapping_quality < MIN_MQ:
                        continue
                    rb = merged_bases([read])
                    if all(rb.get(pos) == base for pos, base in zip(phase_positions, phase_expected)):
                        c['strict_read_329_349_exact'] += 1
    rows, base_rows, summary = [], [], {}
    for name, target in targets.items():
        c = counts[name]
        for i, expected_aa in enumerate(target['protein']):
            pos = target['first']+i
            counts_aa = c['codons'][i]
            sites = target['coding_positions_0based'][3*(pos-1):3*pos]
            row = dict(dataset=kind, element=name, allele=target['allele'],
                       domain_aa=i+1, protein_aa=pos, expected_aa=expected_aa,
                       reference_nt_positions=','.join(str(p+1) for p in sites),
                       expected_aa_fragments=counts_aa[expected_aa],
                       other_aa_fragments=sum(v for aa,v in counts_aa.items() if aa != expected_aa),
                       observed_amino_acids=json.dumps(dict(counts_aa), sort_keys=True),
                       strict_mapping_expected_aa_fragments=c['strict_codons'][i][expected_aa],
                       expected_aa_distinct_fragment_coordinates=len(c['starts'][i][expected_aa]),
                       min_high_quality_base_depth=min(sum(b.values()) for b in c['bases'][3*i:3*i+3]))
            rows.append(row)
        lo, hi = 3*(target['first']-1), 3*target['last']
        positions = target['coding_positions_0based'][lo:hi]
        sequence = target['coding_sequence'][lo:hi]
        for i, (pos, base) in enumerate(zip(positions, sequence)):
            b = c['bases'][i]
            base_rows.append(dict(dataset=kind, element=name, allele=target['allele'],
                                  domain_nt=i+1, reference_nt=pos+1, expected_base=base,
                                  expected_base_fragments=b[base],
                                  alternative_base_fragments=sum(v for x,v in b.items() if x != base),
                                  observed_bases=json.dumps(dict(b),sort_keys=True),
                                  expected_is_unique_majority=bool(b[base] and all(b[base] > v for x,v in b.items() if x != base))))
        group = rows[-len(target['protein']):]
        summary[name] = dict(allele=target['allele'], protein=target['protein'],
                            protein_interval_1based=[target['first'],target['last']],
                            min_expected_codon_fragments=min(r['expected_aa_fragments'] for r in group),
                            max_expected_codon_fragments=max(r['expected_aa_fragments'] for r in group),
                            min_distinct_fragment_coordinates=min(r['expected_aa_distinct_fragment_coordinates'] for r in group),
                            min_strict_mapping_expected_codon_fragments=min(r['strict_mapping_expected_aa_fragments'] for r in group),
                            every_nucleotide_majority_matches_reference=all(r['expected_is_unique_majority'] for r in base_rows[-len(sequence):]),
                            min_expected_nucleotide_fragments=min(r['expected_base_fragments'] for r in base_rows[-len(sequence):]),
                            codons_without_expected_support=[r['protein_aa'] for r in group if not r['expected_aa_fragments']],
                            codons_with_alternative_amino_acids=[r for r in group if r['other_aa_fragments']],
                            qualifying_fragments=c['qualifying_fragments'],
                            full_domain_exact_fragments=c['full_domain_exact_fragments'],
                            full_domain_distinct_coordinates=len(c['distinct_full_domain_fragment_coordinates']),
                            phase_329_349_exact_fragments=c['phase_329_349_exact_fragments'],
                            strict_read_329_349_exact=c['strict_read_329_349_exact'])
    return dict(totals=dict(totals), targets=summary), rows, base_rows


def exact_read_screen(f1, f2, details, all_references):
    """Mapping-independent exact matching of the complete SEC coding stretch.

    Scan forward and reverse-complement reads, with all 78 bases Q30 or higher.
    Also counts a contiguous RNA MITD marker stretch containing its two residues
    that differ from the other called HLA-B allele (positions 329 and 349).
    """
    patterns = {}
    for name, target in details.items():
        first, last = (target['first'], target['last']) if name.startswith('SEC') else (329,349)
        lo, hi = 3*(first-1), 3*last
        seq = target['coding_sequence'][lo:hi]
        patterns[name] = dict(sequence=seq, reverse_complement=str(Seq(seq).reverse_complement()),
                              length_nt=len(seq), amino_acid_interval=[first,last],
                              exact_read_fragments=0, paired_observations=0,
                              competing_reference_names=[r.id for r in all_references if seq in str(r.seq) or str(Seq(seq).reverse_complement()) in str(r.seq)])
    old = patterns['SEC_C0102']['sequence']
    assert old[72:75] == 'TGC'
    old = old[:72] + 'GGC' + old[75:]
    patterns['SEC_current_G25'] = dict(sequence=old,reverse_complement=str(Seq(old).reverse_complement()),
                                     length_nt=len(old),amino_acid_interval=[1,26],
                                     exact_read_fragments=0,paired_observations=0,
                                     competing_reference_names=[r.id for r in all_references if old in str(r.seq)])
    lengths = collections.Counter()
    n = 0
    with f1.open() as a, f2.open() as b:
        while True:
            h1, h2 = a.readline().strip(), b.readline().strip()
            if not h1 and not h2:
                break
            s1, s2 = a.readline().strip(), b.readline().strip()
            a.readline(); b.readline()
            q1, q2 = a.readline().strip(), b.readline().strip()
            assert h1.split()[0].removesuffix('/1') == h2.split()[0].removesuffix('/2'), (h1,h2)
            n += 1
            lengths[len(s1)] += 1; lengths[len(s2)] += 1
            for p in patterns.values():
                hits = 0
                for seq, qual in [(s1,q1),(s2,q2)]:
                    hit = False
                    for needle in [p['sequence'], p['reverse_complement']]:
                        offset = seq.find(needle)
                        while offset >= 0:
                            if min(ord(q)-33 for q in qual[offset:offset+len(needle)]) >= MIN_BQ:
                                hit = True
                            offset = seq.find(needle, offset+1)
                    hits += int(hit)
                p['exact_read_fragments'] += int(hits > 0)
                p['paired_observations'] += hits
    return dict(candidate_fragments=n, read_length_counts=dict(sorted(lengths.items())), patterns=patterns)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reuse-alignments',action='store_true',help='Reuse BAMs from an earlier run of the same reference/input mapping')
    args=parser.parse_args()
    choices = read_families()
    summary = dict(scope='Read evidence only; no construct changes',
                   base_quality_min=MIN_BQ, at_least_one_mate_mapping_quality_min=MIN_MQ,
                   coordinate_deduplication='Reported separately; not UMI molecule counts',
                   reference_policy='One representative per called two-field family, union of DNA and RNA calls; all loci serve as competitors',
                   tools=dict(biopython=Bio.__version__, pysam=pysam.__version__),
                   called_reference_choices=choices, datasets={}, input_sha256={})
    rows, base_rows = [], []
    rna_details, rna_reference = None, None
    for kind, prefix in [('rna','TN25-112103'),('dna','TN25-112103-DNA')]:
        fasta, details, aliases = prepare_reference(kind, choices)
        f1, f2 = [SOURCE / kind / f'{prefix}_candidate_{i}.fq' for i in [1,2]]
        for file in [fasta, f1, f2]:
            summary['input_sha256'][str(file)] = sha(file)
        sam, bam = OUT/f'{kind}.sam', OUT/f'{kind}.name-sorted.bam'
        if not args.reuse_alignments:
            run(['bwa','index',str(fasta)])
            with sam.open('w') as out:
                run(['bwa','mem','-t','8','-a',str(fasta),str(f1),str(f2)],stdout=out)
            run(['samtools','sort','-n','-@','4','-o',str(bam),str(sam)])
            sam.unlink()
        else:
            prior=json.loads((OUT/'read-evidence.json').read_text())
            for file in [fasta,f1,f2]:
                assert prior['input_sha256'][str(file)] == sha(file), 'Cannot reuse alignments for changed inputs'
            assert bam.exists()
        result, values, bases = analyse(kind,bam,details)
        result['reference_aliases'] = aliases
        summary['datasets'][kind] = result
        rows.extend(values)
        base_rows.extend(bases)
        (OUT/f'{kind}-target-reference-details.json').write_text(json.dumps(details,indent=2)+'\n')
        if kind == 'rna':
            rna_details, rna_reference = details, fasta
            summary['exact_RNA_read_screen'] = exact_read_screen(f1,f2,details,list(SeqIO.parse(fasta,'fasta')))
    with (OUT/'codon-evidence.tsv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t')
        writer.writeheader(); writer.writerows(rows)
    with (OUT/'nucleotide-evidence.tsv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(base_rows[0]),delimiter='\t')
        writer.writeheader(); writer.writerows(base_rows)
    (OUT/'read-evidence.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2),flush=True)


if __name__ == '__main__':
    main()
