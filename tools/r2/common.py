"""Read the same verified dataset over public HTTPS or authenticated S3."""
import hashlib,json,os,re,time,sys
from pathlib import Path
from urllib.error import HTTPError
from urllib.parse import quote,urlsplit
from urllib.request import Request,urlopen

ACCESS_CONFIG=Path(__file__).resolve().parents[2]/'data/catalog/access.json'
REPO=ACCESS_CONFIG.parents[2]
sys.path.insert(0,str(REPO/'src'))
from yash_vaccine import catalog

def access_config():
    return json.loads(ACCESS_CONFIG.read_text())

def load_environment_file(path=None):
    explicit=path is not None
    path=Path(path or Path.home()/'.config/yash-r2.env').expanduser()
    if path.exists():
        for line in path.read_text().splitlines():
            if line.strip() and not line.lstrip().startswith('#') and '=' in line:
                key,value=line.split('=',1)
                if key.strip().startswith('R2_'):os.environ.setdefault(key.strip(),value.strip().strip('"').strip("'"))
    elif explicit:
        raise FileNotFoundError('Credential file does not exist: '+str(path))

def load_credentials(path=None):
    load_environment_file(path)
    settings=access_config()
    os.environ.setdefault('R2_BUCKET',settings['bucket']);os.environ.setdefault('R2_ENDPOINT',settings['s3_endpoint'])
    for key in ['R2_ACCESS_KEY_ID','R2_SECRET_ACCESS_KEY','R2_ENDPOINT','R2_BUCKET']:
        if not os.environ.get(key):raise ValueError('Missing credential setting: '+key)
    assert os.environ['R2_BUCKET']=='yash-rna-foundry'

def client():
    import boto3
    from botocore.config import Config
    return boto3.client('s3',endpoint_url=os.environ['R2_ENDPOINT'],region_name='auto',
        aws_access_key_id=os.environ['R2_ACCESS_KEY_ID'],aws_secret_access_key=os.environ['R2_SECRET_ACCESS_KEY'],
        config=Config(signature_version='s3v4',max_pool_connections=24,retries={'mode':'adaptive','max_attempts':10},
                      request_checksum_calculation='when_required',response_checksum_validation='when_required'))

class HTTPBody:
    def __init__(self,response):self.response=response
    def read(self,*args):return self.response.read(*args)
    def close(self):self.response.close()
    def iter_chunks(self,chunk_size):
        while chunk:=self.response.read(chunk_size):yield chunk

class PublicDataClient:
    """Read-only object access. No credentials or authorization headers are sent."""
    transport='https'
    def __init__(self,base_url,bucket):
        parsed=urlsplit(base_url)
        local=parsed.hostname in ['localhost','127.0.0.1','::1']
        if parsed.scheme!='https' and not (parsed.scheme=='http' and local):
            raise ValueError('Public downloads require HTTPS (HTTP allowed only on localhost).')
        if not parsed.netloc or parsed.query or parsed.fragment or parsed.username or parsed.password:
            raise ValueError('Invalid public data base URL')
        self.base_url=base_url.rstrip('/');self.bucket=bucket
    def get_object(self,*,Bucket,Key,Range=None):
        if Bucket!=self.bucket or Key.startswith('/') or '..' in Key.split('/'):
            raise ValueError('Unexpected bucket or object key')
        url=self.base_url+'/'+quote(Key,safe='/')
        headers={'User-Agent':'yash-vaccine-data/1'}
        if Range is not None:
            if not re.fullmatch(r'bytes=(?:\d+-\d*|-\d+)',Range):
                raise ValueError('Expected one HTTP byte range')
            headers['Range']=Range
        for attempt in range(4):
            try:
                response=urlopen(Request(url,headers=headers),timeout=120)
                if Range is not None and response.status!=206:
                    response.close()
                    raise ValueError('Download server did not accept the requested byte range; partial data was not appended.')
                return {'Body':HTTPBody(response),'ContentRange':response.headers.get('Content-Range')}
            except HTTPError as error:
                status=error.code;error.close()
                if (status==429 or 500<=status<=599) and attempt<3:
                    time.sleep(2**attempt);continue
                if status in [401,403,404]:
                    if status==404 and Key.startswith(('references/','software/')):
                        raise PermissionError(f'Third-party resource is not mirrored: {Key}. Obtain it from its provider; see docs/external-resources.md.') from None
                    raise PermissionError(f'Public data unavailable (HTTP {status}): {Key}. If the bucket is private, use authorized S3 credentials; see tools/r2/README.md.') from None
                raise

def download_client(credentials_file=None,public_base_url=None,anonymous=False):
    settings=access_config()
    public_url=public_base_url or os.environ.get('R2_PUBLIC_BASE_URL')
    if public_url or anonymous:
        return PublicDataClient(public_url or settings['public_base_url'],settings['bucket'])
    load_environment_file(credentials_file)
    if os.environ.get('R2_PUBLIC_BASE_URL'):
        return PublicDataClient(os.environ['R2_PUBLIC_BASE_URL'],settings['bucket'])
    if credentials_file or any(os.environ.get(k) for k in ['R2_ACCESS_KEY_ID','R2_SECRET_ACCESS_KEY']):
        load_credentials(credentials_file);return client()
    return PublicDataClient(settings['public_base_url'],settings['bucket'])

def bucket_name(s3):return getattr(s3,'bucket',None) or os.environ['R2_BUCKET']

def bundled_manifest():
    settings=access_config();data=(ACCESS_CONFIG.parent/settings['bundled_manifest']).read_bytes()
    assert hashlib.sha256(data).hexdigest()==settings['bundled_manifest_sha256'],'Bundled catalog checksum mismatch'
    generated=catalog.public_manifest(REPO)
    assert catalog.encoded(generated)==data,'Generated download manifest is stale; run python tools/catalog.py --write'
    return generated

def get_bytes(s3,key):
    body=s3.get_object(Bucket=bucket_name(s3),Key=key)['Body']
    try:return body.read()
    finally:body.close()

def get_json(s3,key):return json.loads(get_bytes(s3,key))

def get_manifest(s3,index_key='INDEX.json'):
    index=get_json(s3,index_key)
    data=get_bytes(s3,index['object_manifest_key'])
    if hashlib.sha256(data).hexdigest()!=index['object_manifest_sha256']:
        raise ValueError('Object manifest checksum mismatch')
    return index,json.loads(data)

def stream_digest(s3,key,output=None):
    h=hashlib.sha256();size=0
    body=s3.get_object(Bucket=bucket_name(s3),Key=key)['Body']
    try:
        for chunk in body.iter_chunks(chunk_size=8*1024*1024):
            h.update(chunk);size+=len(chunk)
            if output is not None:output.write(chunk)
    finally:body.close()
    return size,h.hexdigest()
