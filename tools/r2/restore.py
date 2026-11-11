"""Restore selected private R2 data into a portable /data-equivalent directory.

Defaults to listing the plan. Add --execute to download. Every downloaded object
is checked against the full SHA-256 manifest before files are restored.
"""
import argparse,hashlib,json,os,re,tarfile
from pathlib import Path
from common import download_client,bucket_name,bundled_manifest,get_manifest

def download_verified(s3,item,part):
    """Resume a partial object, hashing the saved prefix and the received suffix."""
    if not re.fullmatch(r'[0-9a-f]{64}',item['sha256']) or not isinstance(item['bytes'],int) or item['bytes']<0:
        raise ValueError('Invalid object checksum/size in manifest')
    offset=part.stat().st_size if part.exists() else 0
    if offset>item['bytes']:raise ValueError('Partial file exceeds expected object size: '+str(part))
    h=hashlib.sha256()
    if offset:
        with part.open('rb') as source:
            while chunk:=source.read(8*1024*1024):h.update(chunk)
    size=offset
    if offset<item['bytes']:
        request=dict(Bucket=bucket_name(s3),Key=item['key'])
        if offset:request['Range']=f'bytes={offset}-'
        response=s3.get_object(**request);body=response['Body']
        try:
            if offset:
                match=re.fullmatch(r'bytes (\d+)-(\d+)/(\d+)',response.get('ContentRange') or '')
                if not match or tuple(map(int,match.groups()))!=(offset,item['bytes']-1,item['bytes']):
                    raise ValueError('Unexpected resume range; partial data was not appended: '+item['key'])
                print('RESUMING',item['key'],'at',offset,'bytes',flush=True)
            with part.open('ab' if offset else 'wb') as output:
                for chunk in body.iter_chunks(chunk_size=8*1024*1024):
                    if size+len(chunk)>item['bytes']:
                        raise ValueError('Download exceeds the expected object size: '+item['key'])
                    output.write(chunk);h.update(chunk);size+=len(chunk)
        finally:body.close()
    elif not part.exists():part.touch()
    if size!=item['bytes']:
        raise ValueError('Incomplete download retained; rerun to resume: '+item['key'])
    if h.hexdigest()!=item['sha256']:
        raise ValueError('Download checksum mismatch; remove the invalid partial cache file before retrying: '+str(part))

def safe_path(root,relative):
    root=Path(root).resolve();target=root/relative
    if not target.resolve().is_relative_to(root):raise ValueError('Unsafe archive path: '+relative)
    return Path(os.path.normpath(target))

def unpack(path,root):
    root=Path(root).resolve()
    observed={};expected=None;links=[]
    with tarfile.open(path,'r:gz') as archive:
        for entry in archive:
            if entry.name.startswith('_manifests/'):
                expected={x['path'].removeprefix('/data/'):x for x in map(json.loads,archive.extractfile(entry))};continue
            if entry.isdir():continue
            dest=safe_path(root,entry.name)
            if entry.issym():links.append((dest,entry.linkname));continue
            if not entry.isfile():raise ValueError('Unexpected archive member: '+entry.name)
            dest.parent.mkdir(parents=True,exist_ok=True);h=hashlib.sha256();size=0
            temp=dest.with_name(dest.name+'.r2-restore-part')
            with archive.extractfile(entry) as source,temp.open('wb') as output:
                while chunk:=source.read(8*1024*1024):h.update(chunk);output.write(chunk);size+=len(chunk)
            observed[entry.name]=dict(bytes=size,sha256=h.hexdigest())
            if dest.exists():
                with dest.open('rb') as f:existing=hashlib.file_digest(f,'sha256').hexdigest()
                if existing!=h.hexdigest():raise FileExistsError('Different existing file: '+str(dest))
                temp.unlink()
            else:temp.replace(dest)
            dest.chmod(entry.mode & 0o777)
    assert expected is not None,'Missing embedded payload manifest'
    assert set(observed)=={k for k,v in expected.items() if v['type']=='file'}
    for name,got in observed.items():assert all(got[k]==expected[name][k] for k in ['bytes','sha256']),name
    assert {str(dest.relative_to(root)) for dest,_ in links}=={k for k,v in expected.items() if v['type']=='symlink'}
    for dest,link in links:
        assert expected[str(dest.relative_to(root))]['target']==link
        if link.startswith('/data/'):target=safe_path(root,link.removeprefix('/data/'))
        else:target=safe_path(root,str(dest.parent.relative_to(root)/link))
        relative=os.path.relpath(target,dest.parent);dest.parent.mkdir(parents=True,exist_ok=True)
        if dest.is_symlink():assert dest.resolve()==target
        elif dest.exists():raise FileExistsError(str(dest))
        else:dest.symlink_to(relative)
    return len(observed),len(links)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--destination',type=Path,required=True);p.add_argument('--category',choices=['raw','references','analysis','software','workflow','all'],default='all')
    p.add_argument('--execute',action='store_true');p.add_argument('--index',default='INDEX.json');p.add_argument('--credentials-file',type=Path)
    p.add_argument('--public-base-url',help='Override the configured public HTTPS endpoint')
    p.add_argument('--anonymous',action='store_true',help='Use public HTTPS even if private credentials are available')
    p.add_argument('--key',action='append',help='Download an exact object key from the catalog; may be repeated')
    args=p.parse_args()
    if not args.execute and args.index=='INDEX.json':manifest=bundled_manifest();s3=None
    else:
        s3=download_client(args.credentials_file,args.public_base_url,args.anonymous)
        _,manifest=get_manifest(s3,args.index)
    if args.key:
        missing=set(args.key)-{x['key'] for x in manifest['objects']}
        if missing:p.error('Unknown object key: '+', '.join(sorted(missing)))
    items=[x for x in manifest['objects'] if args.category=='all' or x['category']==args.category]
    if args.key:items=[x for x in items if x['key'] in args.key]
    if args.execute and not items:
        p.error('No project-hosted objects match. Third-party code, models and references are obtained upstream; see docs/external-resources.md.')
    print(json.dumps(dict(objects=len(items),download_bytes=sum(x['bytes'] for x in items),destination=str(args.destination),execute=args.execute),indent=2))
    if not args.execute:return
    root=args.destination.resolve();root.mkdir(parents=True,exist_ok=True);cache=root/'_downloads';cache.mkdir(exist_ok=True)
    for item in items:
        download=cache/item['sha256'];part=download.with_suffix('.part')
        if not download.exists():
            download_verified(s3,item,part)
            part.replace(download)
        else:
            with download.open('rb') as f:assert hashlib.file_digest(f,'sha256').hexdigest()==item['sha256']
        if item['kind']=='archive':unpack(download,root)
        else:
            relative=item.get('restore_path') or item['source_path'].removeprefix('/data/')
            dest=safe_path(root,relative);dest.parent.mkdir(parents=True,exist_ok=True)
            if dest.exists():
                with dest.open('rb') as f:assert hashlib.file_digest(f,'sha256').hexdigest()==item['sha256'],str(dest)
            else:
                # Same-filesystem link avoids a second copy of large BAMs/indexes.
                os.link(download,dest)
        print('RESTORED',item['key'],flush=True)
    print('RESTORE VERIFIED',len(items),'objects',flush=True)

if __name__=='__main__':
    try:main()
    except (PermissionError,ValueError,FileNotFoundError) as error:raise SystemExit(str(error))
