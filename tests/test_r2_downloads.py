"""Exercise public downloads against a local HTTP server, with no S3 dependency."""
import functools,gzip,hashlib,io,json,os,subprocess,sys,tarfile,tempfile,threading,unittest
from http.server import SimpleHTTPRequestHandler,ThreadingHTTPServer
from pathlib import Path
from unittest.mock import patch

REPO=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(REPO/'tools/r2'))
from common import PublicDataClient,download_client


class PublicDownloads(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name);self.web=self.root/'web';self.web.mkdir();self.authorization=[]
        authorization=self.authorization
        class Handler(SimpleHTTPRequestHandler):
            def do_GET(self):
                authorization.append(self.headers.get('Authorization'));super().do_GET()
            def log_message(self,*args):pass
        handler=functools.partial(Handler,directory=str(self.web))
        self.server=ThreadingHTTPServer(('127.0.0.1',0),handler)
        self.thread=threading.Thread(target=self.server.serve_forever,daemon=True);self.thread.start()
        def close():self.server.shutdown();self.server.server_close();self.thread.join()
        self.addCleanup(close);self.url=f'http://127.0.0.1:{self.server.server_port}'
        self.raw_key='raw/TN25-112103/variant + evidence.vcf';data=b'original variant evidence\n'
        path=self.web/self.raw_key;path.parent.mkdir(parents=True);path.write_bytes(data)
        self.raw=dict(key=self.raw_key,category='raw',kind='file',bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),source_path='/data/inputs/TN25-112103/example.vcf')
        archive_path=self.web/'analysis/context.tar.gz';archive_path.parent.mkdir()
        payload=b'ALTERNATE-PROTEIN-CONTEXT\n';member='results/TN25-112103/context.txt'
        ledger=(json.dumps(dict(path='/data/'+member,type='file',bytes=len(payload),sha256=hashlib.sha256(payload).hexdigest()))+'\n').encode()
        with tarfile.open(archive_path,'w:gz') as archive:
            for name,content in [(member,payload),('_manifests/context.jsonl',ledger)]:
                info=tarfile.TarInfo(name);info.size=len(content);archive.addfile(info,io.BytesIO(content))
        encoded=archive_path.read_bytes()
        self.archive=dict(key='analysis/context.tar.gz',category='analysis',kind='archive',bytes=len(encoded),sha256=hashlib.sha256(encoded).hexdigest())
        self.manifest_path=self.web/'manifest.json';self.manifest_path.write_text(json.dumps(dict(objects=[self.raw,self.archive])))
        (self.web/'INDEX.json').write_text(json.dumps(dict(object_manifest_key='manifest.json',object_manifest_sha256=hashlib.sha256(self.manifest_path.read_bytes()).hexdigest())))

    def run_tool(self,tool,*args):
        env={k:v for k,v in os.environ.items() if not k.startswith('R2_')}
        env['HOME']=str(self.root/'empty-home')
        return subprocess.run([sys.executable,'-S',str(REPO/'tools/r2'/tool),*map(str,args)],capture_output=True,text=True,cwd=REPO,env=env)

    def test_restore_and_verify_anonymously(self):
        dest=self.root/'restored'
        result=self.run_tool('restore.py','--destination',dest,'--public-base-url',self.url,'--execute')
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertEqual((dest/'inputs/TN25-112103/example.vcf').read_bytes(),(self.web/self.raw_key).read_bytes())
        self.assertEqual((dest/'results/TN25-112103/context.txt').read_bytes(),b'ALTERNATE-PROTEIN-CONTEXT\n')
        output=self.root/'verified.json'
        result=self.run_tool('verify_remote.py','--public-base-url',self.url,'--output',output)
        self.assertEqual(result.returncode,0,result.stderr)
        report=json.loads(output.read_text());self.assertEqual(report['full_get_sha256_verified_objects'],2)
        self.assertEqual(report['status'],'VERIFIED');self.assertTrue(self.authorization)
        self.assertTrue(all(x is None for x in self.authorization))

    def test_key_selection_and_corrupt_payload(self):
        dest=self.root/'selected'
        result=self.run_tool('restore.py','--destination',dest,'--key',self.raw_key,'--public-base-url',self.url,'--execute')
        self.assertEqual(result.returncode,0,result.stderr);self.assertFalse((dest/'results').exists())
        (self.web/self.raw_key).write_bytes(b'tampered data')
        bad=self.root/'bad'
        result=self.run_tool('restore.py','--destination',bad,'--key',self.raw_key,'--public-base-url',self.url,'--execute')
        self.assertNotEqual(result.returncode,0);self.assertFalse((bad/'inputs/TN25-112103/example.vcf').exists())

    def test_manifest_checksum_and_unavailable_public_access(self):
        self.manifest_path.write_text('{}')
        result=self.run_tool('restore.py','--destination',self.root/'bad-manifest','--public-base-url',self.url,'--execute')
        self.assertNotEqual(result.returncode,0);self.assertIn('manifest checksum mismatch',result.stderr)
        (self.web/'INDEX.json').unlink()
        result=self.run_tool('restore.py','--destination',self.root/'missing','--public-base-url',self.url,'--execute')
        self.assertNotEqual(result.returncode,0);self.assertIn('Public data unavailable',result.stderr)
        self.assertIn('authorized S3 credentials',result.stderr)

    def test_plan_works_without_network_or_credentials(self):
        result=self.run_tool('restore.py','--destination',self.root/'plan','--category','raw')
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertEqual(json.loads(result.stdout)['objects'],10)
        self.assertFalse((self.root/'plan').exists());self.assertEqual(self.authorization,[])

    def test_explicit_public_access_does_not_load_credentials(self):
        with patch('common.load_environment_file',side_effect=AssertionError('should not read secrets')):
            client=download_client(public_base_url=self.url)
        self.assertIsInstance(client,PublicDataClient)

    def test_default_transport_depends_on_available_credentials(self):
        with patch.dict(os.environ,{},clear=True),patch('common.load_environment_file'):
            self.assertIsInstance(download_client(),PublicDataClient)
        marker=object()
        with patch.dict(os.environ,{'R2_ACCESS_KEY_ID':'fixture-id','R2_SECRET_ACCESS_KEY':'fixture-secret'},clear=True),patch('common.load_environment_file'),patch('common.client',return_value=marker):
            self.assertIs(download_client(),marker)


if __name__=='__main__':unittest.main()
