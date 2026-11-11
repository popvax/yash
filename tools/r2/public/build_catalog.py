"""Prepare the patient/results-only read allowlist and virtual public manifest.

Generation does not deploy the Worker or change the private bucket's permissions.
"""
from pathlib import Path
import argparse,hashlib,json
ROOT=Path(__file__).resolve().parents[3]
HERE=Path(__file__).resolve().parent

def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--execute',action='store_true',help='Write local Worker modules')
 args=parser.parse_args()
 access=json.loads((ROOT/'data/catalog/access.json').read_text())
 blob=(ROOT/'data/catalog'/access['bundled_manifest']).read_bytes()
 assert hashlib.sha256(blob).hexdigest()==access['bundled_manifest_sha256']
 source=json.loads(blob)['objects']
 assert all(o['category'] in ('raw','analysis') for o in source),'Public catalog must not mirror third-party software/models/references'
 allowed={o['key']:dict(bytes=o['bytes'],sha256=o['sha256']) for o in source}
 assert len(allowed)==len(source)==145
 assert not any(key.endswith(('.py','.R','.sh','.sif','.pt','.pth','.pkl','.npz')) for key in allowed),'Source code and pretrained models must come from their providers'
 index=dict(schema_version=1,case='TN25-112103',lead_design='design-2',repository_url='https://github.com/popvax/yash',snapshot_date='2026-11-09',private_bucket=True,
  public_access='patient_and_case_analysis_allowlist',object_manifest_key='public/object-manifest.json',object_manifest_sha256=access['bundled_manifest_sha256'],catalogued_data_objects=len(source),
  original_normal_FASTQs_included=False,third_party_software_models_references_mirrored=False,
  external_resources_url='https://github.com/popvax/yash/blob/main/docs/external-resources.md',
  verification_receipt_url='https://github.com/popvax/yash/blob/main/data/catalog/public-access.json',
  third_party_notices_url='https://github.com/popvax/yash/blob/main/THIRD_PARTY_NOTICES.md')
 print(json.dumps(dict(execute=args.execute,data_objects=len(source),third_party_objects_allowed=0,metadata_objects=0,virtual_objects=2,bucket_permissions_changed=False)))
 if args.execute:
  (HERE/'objects.mjs').write_text('// Generated: patient sequencing and case-specific results only.\nexport default '+json.dumps(allowed,sort_keys=True,indent=2)+';\n')
  (HERE/'index.mjs').write_text('// Generated public entry point; private native INDEX is unchanged.\nexport default '+json.dumps(index,indent=2)+';\n')
  (HERE/'manifest.mjs').write_text('// Generated exact public-manifest bytes, served without changing native R2 metadata.\nexport default '+json.dumps(blob.decode())+';\n')
if __name__=='__main__':main()
