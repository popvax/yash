"""Check selection semantics that would otherwise silently change the antigen pool."""
import csv,json,subprocess,sys,tempfile,unittest
from pathlib import Path
class AdapterTests(unittest.TestCase):
 def test_positive_zero_missing_and_original_measurements(self):
  with tempfile.TemporaryDirectory() as t:
   d=Path(t)
   (d/'expr').write_text('##Gene\tTPM\nPOS\t12.5\nZERO\t0\n')
   (d/'genes').write_text('POS\nZERO\nMISSING\n')
   (d/'quant').write_text('Name\tLength\tEffectiveLength\tTPM\tNumReads\nT1.2\t100\t80\t999\t123\nT2.1\t100\t80\t999\t123\nT3.1\t100\t80\t999\t123\n')
   (d/'gtf').write_text(''.join(f'chr1\tx\ttranscript\t1\t100\t.\t+\t.\tgene_name "{g}"; transcript_id "{tx}"; transcript_type "protein_coding";\n' for g,tx in [('POS','T1.2'),('ZERO','T2.1'),('MISSING','T3.1')]))
   args=[sys.executable,str(Path(__file__).with_name('prepare_expression.py'))]
   for k,v in [('expression','expr'),('quants','quant'),('gtf','gtf'),('genes','genes'),('out','out')]:args += ['--'+k,str(d/v)]
   subprocess.run(args,check=True,capture_output=True)
   with (d/'out/expression-adapted.quant.sf').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
   self.assertEqual(len(rows),1);self.assertEqual(rows[0]['Name'],'T1.2')
   self.assertEqual(rows[0]['TPM'],'999');self.assertEqual(rows[0]['NumReads'],'123')
   self.assertEqual(float(rows[0]['SuppliedGeneExpression']),12.5)
   with (d/'out/cta-gene-selection-audit.tsv').open() as f:audit={r['gene']:r for r in csv.DictReader(f,delimiter='\t')}
   self.assertEqual(audit['ZERO']['status'],'zero');self.assertEqual(audit['MISSING']['status'],'missing')
   self.assertEqual(json.loads((d/'out/expression-adapter-summary.json').read_text())['selected_transcripts'],1)
if __name__=='__main__':unittest.main()
