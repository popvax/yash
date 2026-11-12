"""Native GATK diagnostic: retain global priors but omit segmentation.

This tests dependence on the optional allele-balance model. It does not select
this model, assign somatic status, or replace the recorded segmented workflow.
"""
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import pysam

root=Path(sys.argv[1])
request=json.loads((root/'audit-request.json').read_text())
command=next(c for c in json.loads((root/'original-filter-commands.json').read_text()) if 'FilterMutectCalls' in c)
command=list(command)
at=command.index('--tumor-segmentation')
del command[at:at+2]
command[command.index('-O')+1]=str(root/'diploid-filtered.vcf.gz')
command[command.index('--filtering-stats')+1]=str(root/'diploid-filtering-stats.tsv')
(root/'diploid-command.json').write_text(json.dumps(dict(command=command,
    worker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    change='Omit optional --tumor-segmentation only; native default germline heterozygote expectation is 0.5. Global priors, full input call set, contamination, orientation and FDR unchanged.'),indent=2)+'\n')
with (root/'diploid.stdout.log').open('w') as so,(root/'diploid.stderr.log').open('w') as se:
    subprocess.run(command,check=True,stdout=so,stderr=se)
targets={tuple(r['allele']):r['variant_id'] for r in request['mutation_prior_sensitivity']}
selected=set(request['selected_mutation_ids'])
result=[]
for r in pysam.VariantFile(str(root/'diploid-filtered.vcf.gz')):
    for i,alt in enumerate(r.alts):
        key=(r.chrom,r.pos,r.ref,alt)
        if key in targets:
            result.append(dict(variant_id=targets[key],selected=targets[key] in selected,
                FILTER=list(r.filter),GERMQ=r.info.get('GERMQ'),
                AS_FilterStatus=r.info.get('AS_FilterStatus'),POPAF=r.info.get('POPAF')))
assert selected<=set(r['variant_id'] for r in result)
(root/'diploid-results.json').write_text(json.dumps(result,indent=2)+'\n')
print('DIPLOID_DIAGNOSTIC_COMPLETE',len(result),flush=True)
