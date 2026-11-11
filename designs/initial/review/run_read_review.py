from pathlib import Path
import argparse
import json

p=argparse.ArgumentParser();p.add_argument('--execute',action='store_true');a=p.parse_args()
local=Path(__file__).parent
remote='/data/results/TN25-112103/manual-antigen-review-20261102'
print(json.dumps(dict(execute=a.execute,remote=remote,changes_selection=False)),flush=True)
if not a.execute:raise SystemExit(0)
import modal
v=modal.Volume.from_name('research-lens-workspace')
with v.batch_upload(force=True) as u:
    for name in ['read_review.py','request.json']:u.put_file(str(local/name),remote.removeprefix('/data/')+'/'+name)
s=modal.Sandbox.create('python','-u',remote+'/read_review.py',remote,
    app=modal.App.lookup('yash-antigen-manual-review',create_if_missing=True),
    image=modal.Image.from_id('im-hX1h6HknDMvasm9Z6q530H'),cpu=2,memory=4096,timeout=1200,volumes={'/data':v})
(local/'modal-state.json').write_text(json.dumps(dict(sandbox_id=s.object_id,remote=remote),indent=2)+'\n')
print(s.object_id,flush=True)
stdout=s.stdout.read();stderr=s.stderr.read();s.wait()
(local/'stdout.log').write_text(stdout);(local/'stderr.log').write_text(stderr)
print(stdout,stderr,flush=True)
assert s.poll()==0,s.poll()
(local/'raw-read-review.json').write_bytes(b''.join(v.read_file(remote.removeprefix('/data/')+'/raw-read-review.json')))
