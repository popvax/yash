"""Plan the data launch; --publish explicitly enables public patient-data reads."""
from pathlib import Path
import argparse,json,subprocess,sys,time

ROOT=Path(__file__).resolve().parents[3]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--publish',action='store_true',help='Make the selected R2 allowlist public, verify it and record availability')
    args=parser.parse_args()
    deploy=['npx','--yes','wrangler@4.157.0','deploy','--config','tools/r2/public/wrangler.launch.jsonc']
    verify=[sys.executable,'tools/r2/public/verify_access.py','--output','data/catalog/public-access.json']
    if not args.publish:
        print('PLAN ONLY. Nothing is deployed or made public.')
        print('Public data URL: https://yash-data.samarth-3f0.workers.dev')
        print('At launch: enable the read-only allowlist, verify anonymous access, then record availability.')
        print('Use --publish only after an explicit public-launch decision.')
        return
    access_path=ROOT/'data/catalog/access.json'
    access=json.loads(access_path.read_text())
    if access['public_base_url']!='https://yash-data.samarth-3f0.workers.dev':
        raise ValueError('Review the changed public endpoint before publishing')
    print('PUBLISHING the selected patient-data allowlist. The bucket remains private.',flush=True)
    try:
        subprocess.run(deploy,cwd=ROOT,check=True)
        for attempt in range(3):
            try:
                subprocess.run(verify,cwd=ROOT,check=True)
                break
            except subprocess.CalledProcessError:
                if attempt==2:raise
                print('Verification has not passed; waiting for route propagation before retrying.',flush=True)
                time.sleep(5)
        receipt=json.loads((ROOT/'data/catalog/public-access.json').read_text())
        if not receipt['public_access_enabled'] or receipt['verification_mode']!='ANONYMOUS_HTTPS':
            raise ValueError('External anonymous verification did not pass')
        access['public_access_enabled']=True
        access['public_access']='Enabled after anonymous HTTPS verification on '+receipt['checked_utc']
        access_path.write_text(json.dumps(access,indent=2)+'\n')
    except (Exception,KeyboardInterrupt):
        try:
            withdrawn=subprocess.run(['npx','--yes','wrangler@4.157.0','deploy','--config','tools/r2/public/wrangler.jsonc'],cwd=ROOT)
        except (OSError,KeyboardInterrupt):
            raise SystemExit('Launch failed and withdrawal was not confirmed. Inspect the Worker immediately.')
        if withdrawn.returncode:
            raise SystemExit('Launch verification failed and withdrawal was not confirmed. Inspect the Worker immediately.')
        access['public_access_enabled']=False
        access['public_access']='Disabled after failed launch verification; the disabled Worker profile was restored.'
        access_path.write_text(json.dumps(access,indent=2)+'\n')
        raise SystemExit('Launch verification failed; the disabled Worker profile was restored.')
    print('PUBLIC DATA VERIFIED. Commit the access configuration and dated receipt before making the research repo public.')


if __name__=='__main__':main()
