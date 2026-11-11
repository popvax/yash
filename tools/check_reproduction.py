"""Reproduce lead Design #2 in a disposable directory."""
import json,sys
from pathlib import Path
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'src'))
from yash_vaccine.designs import reproduce

def main():print(json.dumps(reproduce(root),indent=2))

if __name__=='__main__':main()
