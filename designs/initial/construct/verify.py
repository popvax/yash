"""Check the final full-chain scan against independent fresh model predictions."""
import argparse
from pathlib import Path
from common import (HERE, sha, read_json, write_json, read_fasta, rows, inputs,
                    full_scan, load_predictions)
from run import summarize_predictions
from yash_vaccine.hydrophobicity import protein_screen


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--protein", type=Path, default=HERE / "protein")
    p.add_argument("--inference", type=Path, default=HERE / "selected-inference")
    p.add_argument("--output", type=Path, default=HERE / "final-check.json")
    args = p.parse_args()
    if args.output.exists():
        p.error("Choose a new verification output")
    expected = read_json(args.protein / "verification.json")["output_sha256"]
    for name, digest in expected.items():
        if sha(args.protein / name) != digest:
            raise ValueError("Protein export changed: " + name)
    req = read_json(args.protein / "verification-request.json")
    fresh = load_predictions(args.inference, req)
    regions, _, sec, mitd, alleles, _ = inputs()
    sequence = read_fasta(args.protein / "cassette.fasta")
    parts = rows(args.protein / "components.tsv")
    joins = rows(args.protein / "junctions.tsv")
    for join in joins:
        for key in ["junction_index", "linker_start_aa", "linker_end_aa"]:
            join[key] = int(join[key])
        a, b = join["linker_start_aa"], join["linker_end_aa"]
        assert sequence[a:b] == join["linker_sequence"]
        assert sequence[a-10:a] == join["left"][-10:]
        assert sequence[b:b+10] == join["right"][:10]
    assert sequence.startswith(sec) and sequence.endswith(mitd) and len(joins) == len(regions)+1
    cursor = 0
    for part in parts:
        a,b=int(part["start_aa"]),int(part["end_aa"])
        assert a==cursor and sequence[a:b]==part["aa"]
        cursor=b
    assert cursor==len(sequence)
    _,hydrophobicity=protein_screen(sequence,parts)
    assert hydrophobicity['status']=='PASS'
    assert hydrophobicity==read_json(args.protein/'hydrophobicity.json')
    assert sorted((x["region_id"],x["aa"]) for x in parts if x["kind"]=="antigen") == \
           sorted((r["region_id"],r["sequence"]) for r in regions)
    scan = full_scan(sequence, joins)
    assert sorted({r["peptide"] for r in scan}) == req["peptides"]
    cached = rows(args.protein / "boundary-predictions.tsv")
    actual=[]
    affinity_diffs, rank_diffs = [], []
    for row in scan:
        for allele in alleles:
            a,r = fresh[row["peptide"],allele]
            actual.append(dict(**row,allele=allele,affinity_nm=a,affinity_percentile=r))
    assert len(actual)==len(cached)
    for old,new in zip(cached,actual):
        for key in ["junction_index","category","start_aa","end_aa","peptide","allele"]:
            assert str(old[key])==str(new[key]),key
        a,r=float(old["affinity_nm"]),float(old["affinity_percentile"])
        na,nr=new["affinity_nm"],new["affinity_percentile"]
        affinity_diffs.append(abs(na/a-1)); rank_diffs.append(abs(nr-r))
        assert [a<=50,a<=500,r<=.5,r<=2] == [na<=50,na<=500,nr<=.5,nr<=2]
    assert max(affinity_diffs)<1e-4 and max(rank_diffs)<1e-4
    result=dict(status="FRESH_FULL_JUNCTION_CHECK_PASSED",protein_sha256=sha(args.protein/"cassette.fasta"),
        components_sha256=sha(args.protein/"components.tsv"),protein_length_aa=len(sequence),
        all_antigens_exactly_once=True,antigen_regions=len(regions),exact_SEC_and_MITD=True,scored_junctions=len(joins),
        independent_full_protein_scan=True,unique_peptides=len(req["peptides"]),
        unique_peptide_hla_pairs=len(fresh),positional_peptide_hla_rows=len(actual),
        threshold_classifications_identical=True,max_relative_affinity_difference=max(affinity_diffs),
        max_absolute_rank_difference=max(rank_diffs),
        hydrophobicity=hydrophobicity,
        totals=summarize_predictions(actual),by_category=[dict(category=c,**summarize_predictions(
            [r for r in actual if r["category"]==c])) for c in ["SEC","internal","MITD"]],
        input_sha256={"protein_verification":sha(args.protein/"verification.json"),
                      "fresh_inference":sha(args.inference/"inference.json"),
                      "fresh_predictions":sha(args.inference/"predictions.tsv")},
        verification_script_sha256=sha(Path(__file__)),scope="Binding affinity and hydropathy compatibility screen; no experimental topology validation")
    write_json(args.output,result)
    print(__import__("json").dumps(result,indent=2))


if __name__ == "__main__":
    main()
