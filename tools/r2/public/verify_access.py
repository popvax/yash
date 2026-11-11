"""Verify the allowlisted HTTP endpoint without credentials or large downloads."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from urllib.request import Request,urlopen
from urllib.error import HTTPError
from urllib.parse import quote,urlsplit
import argparse,datetime,hashlib,json,sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'tools/r2'))
from common import PublicDataClient,access_config,bundled_manifest,get_manifest


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--public-base-url',help='Override the configured URL; HTTP allowed only on loopback')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();settings=access_config();base=(args.public_base_url or settings['public_base_url']).rstrip('/')
    client=PublicDataClient(base,settings['bucket'])  # Validates URL; never loads secrets.
    expected=bundled_manifest();objects=expected['objects'];bykey={o['key']:o for o in objects}
    index,manifest=get_manifest(client)
    assert index['repository_url']=='https://github.com/popvax/yash','Wrong repository in public index'
    assert index['object_manifest_sha256']==settings['bundled_manifest_sha256'],'Unexpected public manifest'
    assert manifest==expected,'Downloaded manifest differs from the bundled release'
    def request(key,method='GET',headers=None,limit=2_000_000):
        try:
            with urlopen(Request(base+'/'+quote(key,safe='/'),method=method,headers={'User-Agent':'yash-vaccine-data/1',**(headers or {})}),timeout=60) as r:
                data=r.read(limit+1) if method!='HEAD' else b''
                if len(data)>limit:raise ValueError('Response exceeded the bounded verification size: '+key)
                return r.status,r.headers,data
        except HTTPError as error:
            with error:return error.code,error.headers,error.read(min(limit,8192))
    def head(item):
        status,h,_=request(item['key'],'HEAD')
        assert status==200 and int(h['Content-Length'])==item['bytes'] and h['X-Content-SHA256']==item['sha256'],item['key']
    with ThreadPoolExecutor(max_workers=8) as pool:list(pool.map(head,objects))
    native=json.loads((ROOT/'data/catalog/object-manifest.json').read_text())
    denied=[o['key'] for o in native['objects'] if o['key'] not in bykey]
    denied+=['raw/TN25-112103/normal/R1.fastq.gz','.git/config','not-allowed.json']
    def reject(key):assert request(key,'HEAD')[0]==404,'Unlisted key accessible: '+key
    with ThreadPoolExecutor(max_workers=8) as pool:list(pool.map(reject,denied))
    checked=[]
    for key in ['raw/TN25-112103/vendor-transformed-expression.tsv','raw/TN25-112103/DNA_TN25-112103.vcf','analysis/TN25-112103/design-2/paired/filter/pass.vcf.gz']:
        item=bykey[key];assert item['bytes']<2_000_000
        status,h,data=request(key);assert status==200 and len(data)==item['bytes'] and hashlib.sha256(data).hexdigest()==item['sha256'],key
        checked.append(dict(key=key,bytes=len(data),sha256=item['sha256']))
    bam=bykey['raw/TN25-112103/normal/normal.recal.bam'];ranges=[]
    receipt=ROOT/'data/catalog/public-access.json';prior=json.loads(receipt.read_text()) if receipt.exists() else {}
    known=prior.get('normal_bam_range_checks',[])
    for start in [0,bam['bytes']//2,bam['bytes']-65536]:
        end=start+65535;status,h,data=request(bam['key'],headers={'Range':f'bytes={start}-{end}'},limit=65536)
        assert status==206 and len(data)==65536 and h['Content-Range']==f'bytes {start}-{end}/{bam["bytes"]}'
        sha=hashlib.sha256(data).hexdigest();previous=next((r for r in known if r['start']==start),None)
        if previous:assert previous['sha256']==sha,'Range bytes differ from the recorded verification'
        ranges.append(dict(start=start,bytes=len(data),sha256=sha,matched_recorded_range=previous is not None))
    assert request(bam['key'],headers={'Range':f'bytes={bam["bytes"]}-'})[0]==416
    for method in ['POST','PUT','DELETE']:assert request(bam['key'],method)[0]==405
    private=urlsplit(base).hostname in ['localhost','127.0.0.1','::1']
    report=dict(status='VERIFIED',checked_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),public_base_url=base,
                verification_mode='LOOPBACK_PRIVATE_R2_PROXY' if private else 'ANONYMOUS_HTTPS',
                public_access_enabled=not private,credentials_sent=False,bucket_public_access_enabled=False,
                catalogued_data_objects=len(objects),head_checks_passed=len(objects),unlisted_keys_rejected=len(denied),
                full_get_checks=checked,normal_bam_range_checks=ranges,writes_rejected=True,
                full_dataset_redownloaded=False,
                scope='All allowlisted sizes/checksum headers and denied keys; three small full-GET SHA256 checks; three BAM byte ranges. This does not re-download all large objects.')
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(report,indent=2)+'\n')
    print('HTTP ACCESS VERIFIED:',len(objects),'objects;',len(denied),'unlisted keys rejected;',report['verification_mode'])


if __name__=='__main__':
    try:main()
    except (AssertionError,PermissionError,ValueError,FileNotFoundError) as error:raise SystemExit('Verification failed: '+str(error))
