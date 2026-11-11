"""Prepare the recorded analysis image and Ensembl cache; defaults to a plan."""
from pathlib import Path
import argparse,datetime,json
HERE=Path(__file__).resolve().parent
REPO=next(p for p in HERE.parents if (p/'pyproject.toml').is_file())
REMOTE='/data/results/TN25-112103/mutation-evidence-v2-20261102'
PINS=['vaxrank==3.46.0','isovar==1.57.1','varcode==10.16.0','pysam==0.23.3','pyensembl==2.27.3','mhctools==4.0.6','orfdata==0.34.0','gtfparse==3.0.4','datacache==1.27.1']
CODE=r'''
from pathlib import Path
import importlib.metadata,inspect,json,shutil,sys,varcode,isovar
from pyensembl import Genome
from vaxrank.epitope_config import EpitopeConfig
R=Path(sys.argv[1]);R.mkdir(parents=True,exist_ok=True)
verified=R/'analysis-environment-verified.json'
if verified.exists():
 cfg=json.loads((R/'analysis-genome-config.json').read_text());g=Genome(**cfg)
else:
 cfg=json.loads(Path('/data/results/TN25-112103/openvax-mutations-20261018/references/genome-config.json').read_text())
 # datacache uses hard links, unavailable on the Modal volume. Build on local disk.
 cfg['cache_directory_path']='/tmp/yash-v2-ensembl103';g=Genome(**cfg);g.index()
 shutil.copytree(cfg['cache_directory_path'],R/'reference-cache/ensembl103',dirs_exist_ok=True)
 cfg['cache_directory_path']=str(R/'reference-cache/ensembl103')
 (R/'analysis-genome-config.json').write_text(json.dumps(cfg,indent=2)+'\n');g=Genome(**cfg)
v=varcode.Variant('chr3',30688454,'C','CA',genome=g,convert_ucsc_contig_names=True);e=v.effects().top_priority_effect()
assert e.gene_name=='TGFBR2' and type(e).__name__=='FrameShift'
report=dict(versions={n:importlib.metadata.version(n) for n in ['vaxrank','isovar','varcode','pysam','pyensembl','mhctools','mhcflurry','torch','numpy','pandas']},
 annotation_control=dict(gene=e.gene_name,effect=type(e).__name__,short_description=e.short_description),
 isovar_read_collector_signature=str(inspect.signature(isovar.ReadCollector)),protein_sequence_creator_signature=str(inspect.signature(isovar.ProteinSequenceCreator)),epitope_config_signature=str(inspect.signature(EpitopeConfig)))
verified.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2),flush=True)
'''
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--execute',action='store_true');p.add_argument('--output',type=Path,default=REPO/'work/mutation-evidence-v2-20261102');a=p.parse_args()
 print(json.dumps(dict(base_modal_image='im-01M65VX8ESZ8FGAHJD4BHJR6W9',pins=PINS,remote_root=REMOTE,execute=a.execute),indent=2),flush=True)
 if not a.execute:return
 import modal
 a.output.mkdir(parents=True,exist_ok=True)
 image=modal.Image.from_id('im-01M65VX8ESZ8FGAHJD4BHJR6W9').pip_install(*PINS)
 v=modal.Volume.from_name('research-lens-workspace');app=modal.App.lookup('yash-mutation-evidence-v2',create_if_missing=True)
 with modal.enable_output():
  s=modal.Sandbox.create('python','-u','-c',CODE,REMOTE,app=app,image=image,cpu=4,memory=32768,timeout=7200,volumes={'/data':v})
  (a.output/'analysis-image.json').write_text(json.dumps(dict(sandbox_id=s.object_id,image_id=image.object_id,started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),indent=2)+'\n')
  out=s.stdout.read();err=s.stderr.read();s.wait();(a.output/'prepare-environment.log').write_text(out+'\n'+err)
  print(out[-3000:]);print(err[-2000:]);assert s.poll()==0
  (a.output/'analysis-environment-verified.json').write_bytes(b''.join(v.read_file(REMOTE.removeprefix('/data/')+'/analysis-environment-verified.json')))
if __name__=='__main__':main()
