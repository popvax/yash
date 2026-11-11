"""Fresh native class-I affinity ranking after paired DNA and normal-self review."""
import os
os.environ['OMP_NUM_THREADS']='4';os.environ['MKL_NUM_THREADS']='4'
from pathlib import Path
import collections,csv,datetime,hashlib,importlib,importlib.metadata,json,math,sys,time
import numpy as np,pandas as pd,torch
from pyensembl import Genome
from mhctools import MHCflurry
from mhctools.pred import Prediction,Kind
from vaxrank.candidate_epitope import CandidateEpitope
from vaxrank.epitope_config import EpitopeConfig
from vaxrank.epitope_dsl import attach_per_allele_scores,epitopes_for_ranking,default_score_expr
from vaxrank.epitope_logic import slice_epitopes
from vaxrank.vaccine_antigen import VaccineAntigen,TargetableMask,AminoAcidInterval,TumorSpecificityAttestation
from vaxrank.mutant_protein_fragment import MutantProteinFragment
from vaxrank.vaccine_peptide import VaccinePeptide
from vaxrank.reference_proteome import self_reference_matches
from dataclasses import asdict
from self_search import patient_self_matches
R=Path(sys.argv[1]);D=R/'rank';D.mkdir(parents=True,exist_ok=True);torch.set_num_threads(4)
assert json.loads((R/'rna/complete.json').read_text())['status']=='COMPLETED'
proteins=json.loads((R/'rna/primary-proteins.json').read_text());evidence={r['variant_id']:r for r in json.loads((R/'rna/rna-evidence.json').read_text())}
case=json.loads((R/'case.json').read_text());alleles=case['hla_class_I'];assert len(set(alleles))==6
cfg=json.loads((R/'analysis-genome-config.json').read_text());g=Genome(**cfg)
MODELS='/data/results/TN25-112103/openvax-mutations-20261018/references/mhcflurry/2.3.0/models_class1_presentation/models'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def save(name,v):
 p=D/name;t=p.with_suffix('.tmp');t.write_text(json.dumps(v,indent=2,default=str)+'\n');t.replace(p)
def table(name,rows):
 fields=list(dict.fromkeys(k for r in rows for k in r)) or ['variant_id','status']
 with (D/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader()
  for r in rows:w.writerow({k:json.dumps(v) if isinstance(v,(tuple,list,dict)) else v for k,v in r.items()})
def status(phase,**kw):save('status.json',dict(phase=phase,updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),**kw));print(phase,kw,flush=True)
def contribution(nm):
 x=0 if nm>=5000 else (1+math.exp(-350/150))/(1+math.exp((nm-350)/150))
 return x if x>=1e-5 else 0.0
start_time=time.time();ec=EpitopeConfig();expression={};expression_file=Path('/data/inputs/TN25-112103/vendor-transformed-expression.tsv')
model_manifest=json.loads((R/'affinity-model-manifest.json').read_text())
checked_models=[]
for m in model_manifest:
 if not m['path'].startswith('affinity_predictor/'):continue
 path=Path(MODELS)/m['path'];assert path.stat().st_size==m['bytes'] and sha(path)==m['sha256'],str(path)
 checked_models.append(m)
assert checked_models
save('verified-affinity-model-files.json',checked_models)
save('environment-packages.json',{d.metadata['Name']:d.version for d in importlib.metadata.distributions()})
native_modules=['vaxrank.candidate_epitope','vaxrank.epitope_config','vaxrank.epitope_dsl','vaxrank.epitope_logic','vaxrank.vaccine_antigen','vaxrank.mutant_protein_fragment','vaxrank.vaccine_peptide','vaxrank.reference_proteome',MHCflurry.__module__,'mhctools.pred']
save('native-scoring-source-hashes.json',{name:sha(Path(importlib.import_module(name).__file__)) for name in native_modules})
headers=[]
for line in expression_file.read_text().splitlines():
 if line.startswith('##'):headers.append(line);continue
 if not line.strip() or line.startswith('#'):continue
 gene,value=line.split('\t');assert gene not in expression;expression[gene]=float(value)
assert any('WTSPercentileTransformationVersion' in h for h in headers)
save('score-config.json',dict(native_epitope_config=repr(ec),resolved_score_expr=ec.score_expr or default_score_expr(ec),
 peptide_lengths=[8,9,10,11],region_length=25,alleles=alleles,model_bundle='2.3.0',models_path=MODELS,
 weighting='affinity contribution sum * sqrt(supplied gene value * observed RNA alternate fragment fraction)',
 expression_headers=headers,expression_calibration='Vendor-transformed values, not established as ordinary TPM',
 class_II=False,presentation_scoring=False,immunogenicity_prediction=False,selection_quota=None))
occurrences=[]
for p in proteins:
 seq=p['sequence'];lo=p['mutation_start_0based'];hi=p['mutation_end_0based_exclusive']
 assert p['usable_mutant_context'] and 0<=lo<hi<=len(seq) and set(seq)<=set('ACDEFGHIKLMNPQRSTVWY')
 for k in [8,9,10,11]:
  for a in range(len(seq)-k+1):
   if a<hi and a+k>lo:occurrences.append(dict(variant_id=p['variant_id'],protein_id=p['protein_id'],peptide=seq[a:a+k],start=a,end=a+k))
peptides=sorted({o['peptide'] for o in occurrences});table('peptide-occurrences.tsv',occurrences)
affinity={};matches={};normal_matches={};prediction_rows=[]
if peptides:
 status('predicting_classI_affinity',unique_peptides=len(peptides),peptide_HLA_pairs=len(peptides)*len(alleles))
 predictor=MHCflurry(alleles=alleles,default_peptide_lengths=[8,9,10,11],models_path=MODELS)
 raw=predictor.predict_peptides(peptides).to_dataframe();raw.to_csv(D/'native-affinity-predictions.tsv',sep='\t',index=False)
 assert len(raw)==len(peptides)*6 and not raw.duplicated(['peptide','allele']).any()
 assert set(raw.allele)==set(alleles) and np.isfinite(raw.affinity.to_numpy()).all() and (raw.affinity>0).all()
 for r in raw.to_dict('records'):
  affinity[r['peptide'],r['allele']]=float(r['affinity']);prediction_rows.append(dict(peptide=r['peptide'],allele=r['allele'],affinity_nm=float(r['affinity']),percentile_rank=float(r['percentile_rank'])))
 table('peptide-HLA-affinity.tsv',prediction_rows)
 antigen=VaccineAntigen(kind='mutation',amino_acids=proteins[0]['sequence'],
  targetable_mask=TargetableMask((AminoAcidInterval(proteins[0]['mutation_start_0based'],proteins[0]['mutation_end_0based_exclusive']),)),
  tumor_specificity=TumorSpecificityAttestation(status='overridden',evidence_kind='paired_normal_computational_screen',evidence_source='Paired Mutect2, direct normal genotype/read review, population review and Isovar',
   requires_review=True,override_reason='Research prediction with explicit measured normal evidence; not experimentally established tumour specificity'))
 status('exact_reference_peptide_check',peptides=len(peptides))
 matches=self_reference_matches(peptides,antigen,g)
 assert set(matches)==set(peptides) and all(v.source_provenance_complete for v in matches.values())
 save('reference-matches.json',{p:asdict(v) for p,v in matches.items()})
 normal_matches=patient_self_matches(R,peptides)
 save('patient-normal-peptide-matches.json',normal_matches)
occ_by_variant=collections.defaultdict(list)
for o in occurrences:occ_by_variant[o['variant_id']].append(o)
windows=[];ranked=[];exclusions=[];ledger=[];best_hla=[]
for p in proteins:
 vid=p['variant_id'];seq=p['sequence'];lo=p['mutation_start_0based'];hi=p['mutation_end_0based_exclusive'];e=evidence[vid]
 if len(seq)<25:
  exclusions.append(dict(variant_id=vid,gene=p['gene'],reason='RNA_CONTEXT_SHORTER_THAN_25AA'));continue
 eps=[];all_nonself=[]
 for o in occ_by_variant[vid]:
  pep=o['peptide']
  if matches[pep].occurs or normal_matches[pep]['occurs_in_searched_models']:continue
  all_nonself.append(o)
  predictions=tuple(Prediction(kind=Kind.pMHC_affinity,peptide=pep,allele=a,score=affinity[pep,a],value=affinity[pep,a],predictor_name='mhcflurry',predictor_version='2.3.14+models2.3.0') for a in alleles)
  eps.append(CandidateEpitope(sequence=pep,source_sequence=seq,source_name=p['protein_id'],offset=o['start'],predictions=predictions,
   source_class='mutation',overlaps_mutation=True,overlaps_targetable=True,patient_alleles=tuple(alleles),prediction_id=p['protein_id']+':'+str(o['start'])+':'+str(o['end']),
   input_evidence=({'somatic_status':p['somatic_status'],'exact_reference_match':False,'match_in_patient_normal_sequence_models':False},)))
 scored=epitopes_for_ranking(attach_per_allele_scores(eps,ec),ec)
 fragment=MutantProteinFragment(variant=vid,gene_name=p['gene'],amino_acids=seq,mutant_amino_acid_start_offset=lo,mutant_amino_acid_end_offset=hi,
  supporting_reference_transcripts=[],n_overlapping_reads=None,n_alt_reads=None,n_ref_reads=None,n_alt_reads_supporting_protein_sequence=p['supporting_reads'],
  n_overlapping_fragments=e['rna_total_fragments'],n_alt_fragments=e['rna_alt_fragments'],n_ref_fragments=e['rna_ref_fragments'],n_alt_fragments_supporting_protein_sequence=p['supporting_fragments'],
  rna_evidence_method='rna_reads',rna_evidence_subject='fragments')
 antigen=VaccineAntigen(kind='mutation',amino_acids=seq,targetable_mask=TargetableMask((AminoAcidInterval(lo,hi),)),gene_name=p['gene'],transcript_ids=tuple(p['transcript_ids']),source_identifier=vid,
  tumor_specificity=TumorSpecificityAttestation(status='overridden',evidence_kind='paired_normal_computational_screen',evidence_source='Paired DNA and RNA evidence with explicit normal-evidence tiers',requires_review=True,
   override_reason='Computational research ranking with direct normal evidence; uncertainty remains explicit'))
 candidates=[]
 for a,frag in fragment.sorted_subsequences(subsequence_length=25):
  b=a+25
  if not(a<hi and b>lo):continue
  inside=[o for o in all_nonself if a<=o['start'] and o['end']<=b]
  vp=VaccinePeptide(antigen=antigen.sliced(a,b),epitopes=slice_epitopes(scored,a,b),num_target_epitopes_to_keep=0,combined_score_expr='target_epitope_score',ranking_rules=('target_epitope_score',))
  score=float(vp.combined_score);expected=sum(contribution(affinity[o['peptide'],allele]) for o in inside for allele in alleles)
  assert math.isclose(score,expected,rel_tol=1e-8,abs_tol=1e-8),(vid,a,score,expected)
  best=min((affinity[o['peptide'],allele] for o in inside for allele in alleles),default=None)
  value=expression.get(p['gene']);fraction=e['rna_alt_fraction'];proxy=value*fraction if value is not None and value>0 and fraction is not None else None
  entry=dict(variant_id=vid,protein_id=p['protein_id'],gene=p['gene'],DNA_annotation_protein_change=p['DNA_annotation_protein_change'],DNA_annotation_transcript_id=p['DNA_annotation_transcript_id'],RNA_supported_transcript_ids=p['transcript_ids'],sequence=seq[a:b],
   window_start_0based=a,window_end_exclusive=b,target_start_in_window=max(0,lo-a),target_end_in_window=min(25,hi-a),affinity_contribution_sum=score,
   best_target_affinity_nm=best,passes_500nm_gate=best is not None and best<=500,
   supplied_gene_value=value,rna_alternate_fraction=fraction,expression_proxy=proxy,weighted_score=score*math.sqrt(proxy) if proxy is not None else None,
   primary_DNA_eligible=p['primary_DNA_eligible'],gene_assignment_unambiguous=p['gene_assignment_unambiguous'],DNA_evidence_tier=p['DNA_evidence_tier'],somatic_status=p['somatic_status'])
  windows.append(entry);candidates.append(entry)
 if not candidates:
  exclusions.append(dict(variant_id=vid,gene=p['gene'],reason='NO_25AA_WINDOW'));continue
 # Exact maximum after reference exclusion; native centering order breaks ties.
 best=max(candidates,key=lambda r:r['affinity_contribution_sum']);ranked.append(best)
 a=best['window_start_0based'];b=best['window_end_exclusive'];inside=[o for o in all_nonself if a<=o['start'] and o['end']<=b]
 for allele in alleles:
  if not inside:continue
  o=min(inside,key=lambda o:affinity[o['peptide'],allele]);best_hla.append(dict(variant_id=vid,gene=p['gene'],allele=allele,peptide=o['peptide'],affinity_nm=affinity[o['peptide'],allele]))
 for o in inside:
  for allele in alleles:ledger.append(dict(variant_id=vid,peptide=o['peptide'],allele=allele,start=o['start'],end=o['end'],affinity_nm=affinity[o['peptide'],allele],contribution=contribution(affinity[o['peptide'],allele])))
ranked.sort(key=lambda r:(r['weighted_score'] is not None,r['weighted_score'] or 0,r['affinity_contribution_sum']),reverse=True)
primary=[r for r in ranked if r['primary_DNA_eligible'] and r['gene_assignment_unambiguous'] and r['passes_500nm_gate'] and r['weighted_score'] is not None]
for i,r in enumerate(primary,1):r['primary_rank']=i
for i,r in enumerate(ranked,1):r['research_rank']=i
for name,rows in [('all-scored-windows.tsv',windows),('all-ranked-regions.tsv',ranked),('primary-ranked-regions.tsv',primary),('coverage-exclusions.tsv',exclusions),('affinity-contributions.tsv',ledger),('best-peptide-per-HLA.tsv',best_hla)]:table(name,rows)
save('primary-ranked-regions.json',primary)
(D/'primary-regions.fasta').write_text(''.join('>'+r['variant_id']+' '+r['gene']+'\n'+r['sequence']+'\n' for r in primary))
summary=dict(status='COMPLETED',RNA_primary_contexts=len(proteins),unique_peptides=len(peptides),new_peptide_HLA_predictions=len(prediction_rows),exact_reference_matching_peptides=sum(m.occurs for m in matches.values()),
 all_windows=len(windows),ranked_research_regions=len(ranked),primary_ranked_regions=len(primary),primary_genes=len({r['gene'] for r in primary}),
 fixed_quota_applied=False,new_construct_generated=False,normal_DNA_available=True,patient_normal_matching_peptides=sum(m['occurs_in_searched_models'] for m in normal_matches.values()),normal_self_search_proves_normal_absence=False,elapsed_seconds=time.time()-start_time,
 versions={n:importlib.metadata.version(n) for n in ['vaxrank','isovar','varcode','pyensembl','mhctools','mhcflurry','torch']},
 input_sha256=dict(proteins=sha(R/'rna/primary-proteins.json'),evidence=sha(R/'rna/rna-evidence.json'),expression=sha(expression_file),policy=sha(R/'policy.json')),
 outputs={p.name:dict(bytes=p.stat().st_size,sha256=sha(p)) for p in [D/'all-ranked-regions.tsv',D/'primary-ranked-regions.tsv',D/'primary-ranked-regions.json',D/'primary-regions.fasta',D/'best-peptide-per-HLA.tsv']},completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
save('complete.json',summary);status('completed',**summary)
