"""Write a source-backed final comparison after all required stages complete."""
from pathlib import Path
import collections
import csv
import json
import sys

import pysam
import varcode
from pyensembl import Genome
from common import save,receipt,sha
from normal_evidence import NormalEvidence


def read(root,name):return json.loads((root/name).read_text())


def main():
    root=Path(sys.argv[1]);out=root/'report';out.mkdir(parents=True,exist_ok=True)
    stages=['alignment','normal_metadata','bqsr','identity','audit','normal-coding/gather','paired/filter',
            'annotate','rna','normal-self','rank','integrate','construct/protein','construct/verify','construct/mrna']
    receipts={name:read(root,name+'/complete.json') for name in stages}
    assert all(r['status']=='COMPLETED' for r in receipts.values())
    assert receipts['identity']['identity']=='MATCH'
    assert read(root,'construct/final-check.json')['status']=='FRESH_FULL_JUNCTION_CHECK_PASSED'
    old=read(root,'audit/current-mutation-audit.json');by_key=collections.defaultdict(list)
    with pysam.VariantFile(root/'paired/filter/all-unique.vcf.gz') as variants:
        for r in variants:
            by_key[r.contig,r.pos,r.ref,r.alts[0]].append(dict(variant_id=r.id,FILTER=list(r.filter),
                effective_filter=r.info['MN_EFFECTIVE_FILTER'][0],sample_evidence={name:dict(sample.items()) for name,sample in r.samples.items()}))
    annotated={tuple(r['normalized_allele']):r for r in read(root,'annotate/variant-evidence.json')}
    genome=Genome(**read(root,'analysis-genome-config.json'));normal=NormalEvidence(root)
    selected=read(root,'integrate/selected.json');selected_alleles={tuple(r['normalized_allele']) for r in selected if r['source']=='MUTATION'}
    audit=[]
    for row in old:
        key=tuple(row['allele']);calls=by_key.get(key,[]);evidence=annotated.get(key)
        variant=varcode.Variant(*key,genome=genome,convert_ucsc_contig_names=True)
        measured=normal.assess(list(key),variant)
        normal_status=measured['classification']
        if evidence:assert normal_status==evidence['normal_evidence']['classification']
        focal_present=row['normal_verdict']=='TARGET_ALLELE_PRESENT_IN_NORMAL'
        if focal_present and normal_status=='NORMAL_REFERENCE_SUPPORTED':status='NORMAL_GENOTYPE_DISAGREEMENT'
        elif normal_status=='PRESENT_IN_NORMAL' or focal_present:status='PRESENT_IN_NORMAL'
        elif normal_status=='UNRESOLVED_NORMAL_EVIDENCE':status='UNRESOLVED_NORMAL_EVIDENCE'
        elif evidence and evidence['primary_DNA_eligible']:status='SUPPORTED_PAIRED_SOMATIC'
        elif calls and all(r['effective_filter']!='PASS' for r in calls):status='NOT_SUPPORTED_BY_PAIRED_NATIVE_FILTERS'
        elif evidence:status='REQUIRES_NORMAL_POPULATION_OR_CONTEXT_REVIEW'
        elif not calls:status='NORMAL_REFERENCE_SUPPORTED_NOT_REDISCOVERED_AS_EXACT_ALLELE'
        else:status='UNRESOLVED'
        audit.append(dict(row,prior_audit_somatic_classification=row['somatic_classification'],
            somatic_classification=status,final_reassessment=status,paired_exact_allele_records=calls,
            new_normal_evidence=measured,
            reassessed_normal_REF_fragment_count=measured['read_evidence']['fragment_counts'].get('REF',0),
            reassessed_normal_ALT_fragment_count=measured['read_evidence']['fragment_counts'].get('ALT',0),
            new_primary_DNA_eligible=evidence['primary_DNA_eligible'] if evidence else False,
            retained_in_new_selection=key in selected_alleles,
            no_exact_call_does_not_mean_native_failure=True))
    save(out/'current-17-final-audit.json',audit)
    with (out/'current-17-final-audit.tsv').open('w') as stream:
        fields=['antigen_id','gene','protein_change','allele','reassessed_normal_REF_fragment_count','reassessed_normal_ALT_fragment_count','normal_verdict','final_reassessment','retained_in_new_selection']
        writer=csv.DictWriter(stream,fieldnames=fields,delimiter='\t',extrasaction='ignore');writer.writeheader()
        for row in audit:writer.writerow({k:json.dumps(v) if isinstance(v,(list,dict)) else v for k,v in row.items()})
    callable_summary=read(root,'normal-coding/gather/callability-summary.json')
    genotype_scope=read(root,'normal-self/scope.json')
    summary=dict(status='COMPLETED',identity='MATCH',current_mutation_reassessment=dict(collections.Counter(r['final_reassessment'] for r in audit)),
        normal_callability=callable_summary,normal_self_search_scope=genotype_scope,
        stage_receipts=receipts,new_selected_sources=dict(collections.Counter(r['source'] for r in selected)),
        selected_regions=len(selected),new_mutation_genes=sorted({r['gene'] for r in selected if r['source']=='MUTATION'}),
        original_initial_design_modified=False,website_modified=False,pushed_to_public_repository=False)
    save(out/'summary.json',summary)
    normal_pass=callable_summary['genotype_covered_bases']['GQ20_DP20'];target_bases=callable_summary['target_bases']
    lines=['# Matched-normal reassessment','',
        'This analysis uses the new white-blood-cell DNA with the original tumour DNA and RNA. The recorded initial design remains unchanged for comparison.','',
        f"The normal and tumour fingerprints matched. Within the explicitly requested normal coding-plus-flank intervals, {normal_pass:,} of {target_bases:,} bases ({normal_pass/target_bases:.1%}) had native genotype quality at least 20 and depth at least 20. This is not a whole-genome callability percentage.",'',
        '## The 17 previously selected mutations','',
        '| Gene | Protein change | Normal REF / ALT fragments | Reassessment | Selected now |','| --- | --- | --- | --- | --- |']
    for row in audit:lines.append(f"| {row['gene']} | {row['protein_change']} | {row['reassessed_normal_REF_fragment_count']} / {row['reassessed_normal_ALT_fragment_count']} | {row['final_reassessment']} | {row['retained_in_new_selection']} |")
    lines+=['','Missing normal coverage is unresolved. A missing exact allele in independent discovery is not labelled a native filtering failure. Blood-cell presence is not automatically labelled inherited.','',
        '## Independent discovery and RNA reconstruction','',
        'New paired Mutect2 discovery covered the planned primary nuclear chromosome intervals without using Caris calls or the old 17 mutations as a nomination list. Native sample-wide filtering used the matched contamination/segmentation estimates, read-orientation evidence and a 5% model-FDR target; this is not an empirically measured false-discovery rate.',
        '',f"Annotation retained {receipts['annotate']['primary_DNA_eligible']} primary DNA candidates. Isovar reconstructed {receipts['rna']['usable_primary_proteins']} usable leading protein contexts; {receipts['rank']['primary_ranked_regions']} mutation regions passed the explicit primary RNA/expression/class-I ranking gates before construct hydropathy selection.",'',
        '## Selected design and checks','',f"The new selection contains {len(selected)} regions: "+', '.join(f'{count} {source}' for source,count in summary['new_selected_sources'].items())+'.',
        'All passing primary mutation regions were considered; the mutation branch has no fixed top-N quota. The recorded RNA-source discovery and normal-tissue splice screen were reused with their provenance. Altered-target scoring also excludes possible matches in the new normal coding/HLA search models. CTA regions remain intentional unmutated self targets.',
        '', 'Every internal and terminal linker/order option was freshly predicted with the recorded class-I affinity model. The joint optimizer, independent full-chain junction scan, separate fresh selected-junction predictions and hydropathy checks passed. The fixed protein was re-encoded with the recorded PopRNA settings; translation and GenBank annotations were verified. PolyA and the downstream extension remain unspecified, so the exported nucleotide segments are partial.',
        '', '## Practical limits','',
        'Supported paired somatic status is evidence from these samples and models, not a guarantee of healthy-tissue absence. Normal coding models include possible local unphased combinations and explicitly retained unresolved cases; no search establishes a complete personal proteome. Normal DNA does not measure normal RNA splicing, protein translation or HLA presentation. Affinity predictions and synonymous encoding heuristics do not establish immunogenicity, efficacy, topology, expression or manufacturing suitability.',
        '', 'The normal sequencing instrument is unconfirmed. BaseRecalibrator uses its native --default-platform UNKNOWN setting, preserving the original alignment header. The original Caris reference with its custom masks is unavailable; the normal was aligned to the same compatible GRCh38 analysis reference used for this reassessment.',
        '', 'See summary.json for native receipts and callability scope, current-17-final-audit.json for full focal evidence, and the integrate comparison/selection tables for every added, removed and retained region. No website or public repository was changed by this analysis.','']
    (out/'report.md').write_text('\n'.join(lines))
    receipt(out,[out/name for name in ['summary.json','current-17-final-audit.json','current-17-final-audit.tsv','report.md']],worker_sha256=sha(__file__),selected_regions=len(selected),all_required_stage_receipts_passed=True)


if __name__=='__main__':main()
