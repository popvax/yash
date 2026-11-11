"""Prepare a patient/results-only public catalog; third-party resources stay upstream.

Native manifests and private R2 objects are preserved. This only prepares local
publication files; build_catalog.py and an explicit deployment publish the change.
"""
from pathlib import Path
import argparse,hashlib,json
ROOT=Path(__file__).resolve().parents[3]
CAT=ROOT/'data/catalog'

def encoded(x):return (json.dumps(x,indent=2)+'\n').encode()
def digest(b):return hashlib.sha256(b).hexdigest()

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--execute',action='store_true');args=p.parse_args()
 body=(CAT/'object-manifest.json').read_bytes();native=json.loads(body)
 access=json.loads((CAT/'access.json').read_text())
 expected=access.get('native_manifest_sha256',access['bundled_manifest_sha256'])
 if digest(body)!=expected:raise ValueError('Historical object manifest changed')
 def public(x):return x['category'] in ('raw','analysis') and not x['key'].endswith(('.py','.R','.sh','.sif','.pt','.pth','.pkl','.npz'))
 hosted=[x for x in native['objects'] if public(x)]
 external=[x for x in native['objects'] if not public(x)]
 assert len(hosted)==145 and len(external)==315
 public=dict(schema_version=1,case=native['case'],snapshot_date='2026-11-09',
  scope='Patient sequencing and case-specific analysis outputs only. Code lives in GitHub. Third-party software, trained weights and general reference libraries are upstream-only. Normal FASTQs excluded.',
  native_manifest_sha256=digest(body),objects=hosted)
 public_body=encoded(public)
 # Provider pointers are curated acquisition metadata, independent of licensing audits.
 # New or changed resources need explicit provider metadata rather than inferred terms.
 external_doc=json.loads((CAT/'external-resources.json').read_text())
 pointers={x['former_project_key']:x for x in external_doc['objects']}
 if len(pointers)!=len(external) or set(pointers)!={x['key'] for x in external}:
  raise ValueError('Review external-resources.json for the changed native inventory')
 for item in external:
  pointer=pointers[item['key']]
  if any(pointer[k]!=item[k] for k in ('bytes','sha256')) or not pointer.get('provider_urls'):
   raise ValueError('Provider identity missing or changed: '+item['key'])
 access.update(bundled_manifest='public-object-manifest.json',bundled_manifest_sha256=digest(public_body),native_manifest='object-manifest.json',native_manifest_sha256=digest(body))
 print(json.dumps(dict(execute=args.execute,public_objects=len(hosted),upstream_only_objects=len(external),withdrawn_bytes=sum(x['bytes'] for x in external),bucket_objects_deleted=False)))
 if args.execute:
  (CAT/'public-object-manifest.json').write_bytes(public_body)
  (CAT/'external-resources.json').write_bytes(encoded(external_doc))
  (CAT/'access.json').write_bytes(encoded(access))
if __name__=='__main__':main()
