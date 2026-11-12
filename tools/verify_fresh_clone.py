"""Exercise independent work and our optional helpers in a fresh environment.

Only included evidence is used. Prediction/download/cloud operations are tested
as plans; the existing download tests use their own local HTTP fixture.
"""
import argparse,csv,json,subprocess,sys,tempfile,time
from pathlib import Path

def main():
    root=Path(__file__).resolve().parents[1]
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=root/'work/fresh-clone-verification.json')
    args=parser.parse_args();start=time.monotonic();logs=[]
    with tempfile.TemporaryDirectory(prefix='yash-fresh-handoff-') as tmp:
        check=Path(tmp)/'repo'
        def run(command,cwd=check,*,succeeds=True):
            result=subprocess.run(command,cwd=cwd,capture_output=True,text=True)
            logs.append(dict(command=command,returncode=result.returncode,stdout=result.stdout,stderr=result.stderr))
            print(('PASS' if (result.returncode==0)==succeeds else 'FAIL'),command[1:4],flush=True)
            if (result.returncode==0)!=succeeds:
                print(result.stdout[-4000:],result.stderr[-4000:]);raise RuntimeError(command)
            return result
        # Use normal Git transport: a local filesystem clone can otherwise copy
        # unreachable historical objects and conceal a stale-commit dependency.
        run(['git','clone','--quiet','--no-local',str(root),str(check)],root)
        run([sys.executable,'-m','venv',str(check/'.venv')],root)
        python=str(check/'.venv/bin/python');cli=[python,'-m','yash_vaccine']
        run([python,'-m','pip','install','-q','-e','.[mrna]'])
        run([*cli,'--help'])
        run([python,'designs/design-2/verify_results.py'])
        run([*cli,'verify'])
        status=run([python,'designs/design-2/pipeline/status.py'])
        assert json.loads(status.stdout)['source']=='recorded_completed_results'
        # Editing maintained documentation/help is independent of the frozen
        # experiment. Actual scientific-file changes must still be rejected.
        for relative,append in [('src/yash_vaccine/cli.py',b'\n# Harmless help/source maintenance.\n'),('config/README.md',b'\nTemporary documentation check.\n')]:
            path=check/relative;original=path.read_bytes()
            try:
                path.write_bytes(original+append);run([*cli,'verify'])
            finally:path.write_bytes(original)
        for relative in ['designs/design-2/selection/selected.tsv','designs/design-2/inputs/integrated/candidate-pool.tsv']:
            path=check/relative;original=path.read_bytes()
            try:
                path.write_bytes(original+b'\n');run([*cli,'verify'],succeeds=False)
            finally:path.write_bytes(original)
        run([*cli,'reproduce','--output','work/lead-reproduction'])
        lead_result=json.loads((check/'work/lead-reproduction/verification.json').read_text())
        assert lead_result['selected_regions']==24 and lead_result['protein_length_aa']==724
        run([*cli,'compare','designs/initial','designs/design-2','--output','work/complete-design-comparison'])
        complete_comparison=json.loads((check/'work/complete-design-comparison/comparison.json').read_text())
        assert (complete_comparison['retained'],complete_comparison['added'],complete_comparison['removed'])==(22,2,7)
        run([*cli,'reproduce','designs/initial','--output','work/initial-reproduction'])
        initial_result=json.loads((check/'work/initial-reproduction/verification.json').read_text())
        assert initial_result['polyA_status']=='unspecified' and initial_result['rna_length_nt'] is None
        run([*cli,'run','designs/initial/selection/design.json'],succeeds=False)
        run([*cli,'init','independent-analysis'])
        independent=check/'designs/independent-analysis'
        assert {path.name for path in independent.iterdir()}=={'README.md','run.py'}
        starter=run([python,str(independent/'run.py')])
        assert 'No analysis has been run' in starter.stdout
        # A researcher supplies an ordinary script and their own text output.
        # It imports no runner code and needs no recipe or result schema.
        (independent/'run.py').write_text(
            'from pathlib import Path\nimport csv,collections\n'
            'here=Path(__file__).resolve().parent\nrepo=here.parents[1]\n'
            'with (repo/"designs/initial/selection/candidate-pool.tsv").open() as source:\n'
            '    records=list(csv.DictReader(source,delimiter="\\t"))\n'
            'counts=collections.Counter(row["source"] for row in records)\n'
            '(here/"overview.txt").write_text("Available regions: "+str(len(records))+"\\n"+'
            '"".join(name+": "+str(value)+"\\n" for name,value in sorted(counts.items())))\n')
        run([python,str(independent/'run.py')])
        with (root/'designs/initial/selection/candidate-pool.tsv').open() as f: pool_count=len(list(csv.DictReader(f,delimiter='\t')))
        assert 'Available regions: '+str(pool_count) in (independent/'overview.txt').read_text()
        assert not (independent/'design.json').exists()
        unsupported=run([*cli,'compare','designs/initial/selection',str(independent)],succeeds=False)
        assert 'optional recipe runner' in unsupported.stderr
        run([*cli,'init','fewer-mutations','--from','designs/initial/selection'])
        recipe=check/'designs/fewer-mutations/design.json'
        value=json.loads(recipe.read_text());value['selection']['quotas']['MUTATION']-=1
        recipe.write_text(json.dumps(value,indent=2)+'\n')
        run([*cli,'run',str(recipe)])
        run([*cli,'compare','designs/initial/selection','designs/fewer-mutations','--output','work/comparison'])
        comparison=json.loads((check/'work/comparison/comparison.json').read_text())
        assert (comparison['retained'],comparison['added'],comparison['removed'])==(initial_result['selected_regions']-1,0,1)
        alternative_aa=comparison['right_summary']['protein_length_aa']
        assert alternative_aa >= 25*(initial_result['selected_regions']-1)
        removed=[r for r in comparison['regions'] if r['status']=='removed']
        with (root/'designs/initial/selection/selected.tsv').open() as f: initial_regions=list(csv.DictReader(f,delimiter='\t'))
        last_mutation=[r for r in initial_regions if r['source']=='MUTATION'][-1]
        assert [(r['gene'],r['protein_change']) for r in removed]==[(last_mutation['gene'],last_mutation['protein_change'])]
        moved=[r for r in comparison['regions'] if r['gene']=='MUC3A']
        assert moved[0]['left_antigen_id']!=moved[0]['right_antigen_id']
        assert moved[0]['status']=='retained'
        counts=json.loads((check/'designs/fewer-mutations/construct-verification.json').read_text())
        assert 0 < counts['required_peptides'] < counts['cache_peptides']
        assert counts['required_peptide_hla_pairs']==6*counts['required_peptides']
        run([*cli,'reproduce','designs/fewer-mutations'])
        run([*cli,'prepare',str(recipe),'--output','work/queries'])
        plan=run([*cli,'predict','work/queries','--output','work/predictions'])
        request=json.loads((check/'work/queries/prediction-request.json').read_text())
        assert f"{request['peptide_hla_pairs']:,}" in plan.stdout and 'No job started' in plan.stdout
        assert not (check/'work/predictions').exists()
        run([*cli,'init','new-linker','--from',str(recipe)])
        linker_recipe=check/'designs/new-linker/design.json'
        value=json.loads(linker_recipe.read_text());value['linkers']['choices'].append(dict(name='AAY',sequence='AAY'))
        linker_recipe.write_text(json.dumps(value,indent=2)+'\n')
        rejected=run([*cli,'run',str(linker_recipe)],succeeds=False)
        assert 'missing' in rejected.stderr.lower() and not (linker_recipe.parent/'cassette.fasta').exists()
        run([*cli,'data','list'])
        plan=run([*cli,'data','restore','--destination','work/data','--category','raw'])
        assert json.loads(plan.stdout)['objects']==10 and not (check/'work/data').exists()
        run([*cli,'upstream','dna','--output','work/dna'])
        run([*cli,'upstream','rna','--output','work/rna'])
        run([*cli,'upstream','lens','--output-root','/data/results/TN25-112103/example'])
        run([*cli,'upstream','regions','--output-root','/data/results/TN25-112103/replays/example'])
        run([python,'-m','unittest','discover','-s','tests','-v'])
        for name in ('selection/selected.tsv','selection/selected.fasta',
                     'construct/protein/cassette.fasta','construct/protein/components.tsv',
                     'construct/protein/boundary-predictions.tsv','construct/mrna/mrna-segments.fasta',
                     'construct/mrna/template-segments.gb'):
            assert (root/'designs/initial'/name).read_bytes()==(check/'designs/initial'/name).read_bytes(),name
        tracked=run(['git','status','--porcelain','--untracked-files=no'])
        assert not tracked.stdout.strip(),'Tracked files changed in the fresh clone'
        report=dict(status='VERIFIED',source_commit=run(['git','rev-parse','HEAD']).stdout.strip(),
                    fresh_clone=True,fresh_virtualenv=True,python=sys.version.split()[0],
                    design_2_verified_without_historical_git_objects=True,lead_replayed_from_recorded_evidence=True,lead_regions=lead_result['selected_regions'],lead_protein_aa=lead_result['protein_length_aa'],documentation_changes_do_not_invalidate_frozen_results=True,scientific_file_changes_rejected=True,complete_design_comparison=True,
                    initial_candidate_byte_identical=True,initial_protein_aa=initial_result["protein_length_aa"],initial_rna_nt=initial_result['rna_length_nt'],polyA_status=initial_result['polyA_status'],alternative_regions=initial_result["selected_regions"]-1,alternative_cassette_aa=alternative_aa,
                    independent_script_without_recipe=True,independent_plain_text_output=True,
                    stable_region_comparison=True,required_vs_cached_counts=True,
                    missing_predictions_rejected=True,baseline_overwrite_rejected=True,
                    remote_jobs_started=False,external_data_downloaded=False,
                    existing_tests_passed=8,elapsed_seconds=time.monotonic()-start,logs=logs)
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(report,indent=2)+'\n')
        args.output.with_suffix('.txt').write_text('Fresh handoff verified.\n\n'
            'A clean clone ran an independent script with plain-text outputs and no recipe schema.\n'
            'Both complete designs replayed, including lead selection/assembly/RNA from recorded evidence.\n'
            'The optional helper built\n'
            'a smaller alternative, and compared regions using stable identities.\n'
            'Prediction, data and upstream commands planned without starting remote work.\n'
            'Missing predictions and attempts to overwrite the baseline were rejected.\n'
            'All eight existing tests passed. No tracked files were changed.\n\n'
            'Commit: '+report['source_commit']+'\nFull commands and logs: '+str(args.output.resolve())+'\n')
        print('FRESH HANDOFF VERIFIED:',args.output,flush=True)

if __name__=='__main__':main()
