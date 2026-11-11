"""Replay the lead design offline from included evidence; no cloud jobs."""
from pathlib import Path
import argparse,json
from yash_vaccine.lead import reproduce
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,help='New directory for replayed outputs')
    args=parser.parse_args()
    print(json.dumps(reproduce(Path(__file__).resolve().parents[2],args.output),indent=2))
