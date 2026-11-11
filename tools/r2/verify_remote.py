"""Independently GET every manifest object and verify its full SHA-256 digest.

Works over public HTTPS or authenticated S3, without the original Modal volume.
Streams into memory-sized buffers, never downloads the dataset to disk.
"""
import argparse,concurrent.futures,json,threading,time
from pathlib import Path
from common import download_client,bucket_name,get_manifest,stream_digest

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--index',default='INDEX.json');p.add_argument('--output',type=Path,required=True)
    p.add_argument('--workers',type=int,default=8);p.add_argument('--credentials-file',type=Path)
    p.add_argument('--resume',type=Path,action='append',default=[],help='Reuse a completed full-GET verification report when key, size and SHA256 all match')
    p.add_argument('--public-base-url',help='Override the configured public HTTPS endpoint')
    p.add_argument('--anonymous',action='store_true',help='Use public HTTPS even if private credentials are available')
    args=p.parse_args();s3=download_client(args.credentials_file,args.public_base_url,args.anonymous);index,manifest=get_manifest(s3,args.index)
    lock=threading.Lock();receipts=[];failures=[];start=time.time();args.output.parent.mkdir(parents=True,exist_ok=True)
    previous={}
    for path in args.resume:
        report=json.loads(path.read_text())
        assert report['status']=='VERIFIED' and not report['failures'] and report['bucket']==bucket_name(s3),str(path)
        assert report['full_get_sha256_verified_objects']==len(report['objects']),str(path)
        for record in report['objects']:
            assert record['status']=='VERIFIED_BY_FULL_GET',record['key']
            previous[record['key']]=record
    pending=[]
    for item in manifest['objects']:
        old=previous.get(item['key'])
        if old and all(old[k]==item[k] for k in ['key','bytes','sha256']):receipts.append(dict(old))
        else:pending.append(item)
    reused=len(receipts)
    def verify(item):
        size,sha=stream_digest(s3,item['key'])
        assert size==item['bytes'] and sha==item['sha256'],item['key']
        with lock:
            receipts.append(dict(key=item['key'],bytes=size,sha256=sha,status='VERIFIED_BY_FULL_GET'))
            status=dict(status='RUNNING',verified_objects=len(receipts),planned_objects=len(manifest['objects']),
                        verified_bytes=sum(x['bytes'] for x in receipts),elapsed_seconds=time.time()-start)
            args.output.with_suffix('.progress.json').write_text(json.dumps(status,indent=2)+'\n')
            print(json.dumps(dict(verified=item['key'],bytes=size,completed=len(receipts),planned=len(manifest['objects']))),flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures={pool.submit(verify,item):item for item in pending}
        for future in concurrent.futures.as_completed(futures):
            try:future.result()
            except Exception as e:failures.append(dict(key=futures[future]['key'],type=type(e).__name__,message=str(e)[:1000]))
    report=dict(status='VERIFIED' if not failures else 'FAILED',bucket=bucket_name(s3),
                manifest_key=index['object_manifest_key'],manifest_sha256=index['object_manifest_sha256'],
                full_get_sha256_verified_objects=len(receipts),planned_objects=len(manifest['objects']),
                verified_bytes=sum(x['bytes'] for x in receipts),elapsed_seconds=time.time()-start,
                verification_inputs='R2 GET responses and SHA256 manifest only',
                reused_full_get_objects=reused,new_full_get_objects=len(receipts)-reused,
                resumed_reports=[str(x) for x in args.resume],failures=failures,objects=sorted(receipts,key=lambda x:x['key']))
    args.output.write_text(json.dumps(report,indent=2)+'\n');print('FINAL',json.dumps({k:v for k,v in report.items() if k!='objects'}),flush=True)
    if failures:raise SystemExit(1)

if __name__=='__main__':
    try:main()
    except (PermissionError,ValueError,FileNotFoundError) as error:raise SystemExit(str(error))
