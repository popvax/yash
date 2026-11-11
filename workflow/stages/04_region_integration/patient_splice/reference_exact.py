"""Complete exact matching against a verified native reference inventory.

The frozen inventory is output from Vaxrank's own protein-source factory.
Validate every association against the current Ensembl GTF/protein FASTA, then
compare a multi-pattern scan with exhaustive k-mer enumeration for every protein.
This changes execution cost, not the exact-substring or source-exclusion policy.
"""
import ast
from collections import defaultdict
import dataclasses
import gzip
import hashlib
import inspect
import json
from pathlib import Path
import re
import sys

import ahocorasick
from Bio import SeqIO
import vaxrank.reference_proteome as native
from vaxrank.vaccine_antigen import SelfReferenceMatch, SelfReferenceSource
from splice_rna import target_peptides

SNAPSHOT = Path("/data/results/TN25-112103/combined-antigen-regions-20261030/native-reference-snapshot.jsonl")
SNAPSHOT_SHA = "25a2f4cbe3baf7f72beeadbcdeafb9288a3c7c1525001cb49098bfe7632d8d23"
BASE = Path("/data/raft/references/homo_sapiens/pyensembl/GRCh38/ensembl103")
PROTEINS = BASE / "Homo_sapiens.GRCh38.pep.all.fa.gz"
GTF = BASE / "Homo_sapiens.GRCh38.103.gtf"
PROTEIN_SHA = "f8df3d7cfc85718e9161579bd9168a3598547ee50022b5e1c272ed66293f0505"


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def ast_digest(function):
    node = ast.parse(inspect.getsource(function)).body[0]
    return hashlib.sha256(ast.dump(node, include_attributes=False).encode()).hexdigest()


def main(root):
    assert digest(SNAPSHOT) == SNAPSHOT_SHA
    assert digest(PROTEINS) == PROTEIN_SHA
    assert ast_digest(native._protein_source_snapshot) == "f2435753c5274f11938313b87bdcb1e24c07ba6f378f2c0e6624f26d6c7acb1c"
    assert ast_digest(native.self_reference_matches) == "864578418c97c26449581fa8d589a83f664d7091ca8117f1a1ae67c832a5cbcb"
    with SNAPSHOT.open() as handle:
        inventory = [json.loads(line) for line in handle]
    with gzip.open(PROTEINS, "rt") as handle:
        proteins = {record.id.split(".")[0]: str(record.seq) for record in SeqIO.parse(handle, "fasta")}
    expected = {}
    with GTF.open() as handle:
        for line in handle:
            columns = line.rstrip().split("\t")
            if len(columns) != 9 or columns[2] != "CDS":
                continue
            attrs = dict(re.findall(r'(\w+) "([^\"]*)"', columns[8]))
            protein = attrs.get("protein_id", "").split(".")[0]
            if protein not in proteins:
                continue
            key = (native.normalize_ensembl_gene_id(attrs["gene_id"]),
                   attrs["transcript_id"].split(".")[0], protein)
            value = (attrs.get("gene_name", ""), attrs.get("transcript_biotype", ""))
            assert key not in expected or expected[key] == value
            expected[key] = value
    associations = {}
    for record in inventory:
        for source in record["sources"]:
            key = source["gene_id"], source["transcript_id"], source["protein_id"]
            assert source["protein_id"] in proteins and record["sequence"] == proteins[source["protein_id"]]
            assert key not in associations
            associations[key] = (source["gene_name"], source["transcript_biotype"])
    assert associations == expected, dict(snapshot=len(associations), current_GTF=len(expected),
          missing=len(set(expected) - set(associations)), extra=len(set(associations) - set(expected)))
    print("REFERENCE_INVENTORY_VERIFIED", len(inventory), "proteins", len(associations), "sources", flush=True)
    data = root / "reconstruction/patient-windows.jsonl.gz"
    with gzip.open(data, "rt") as handle:
        windows = [json.loads(line) for line in handle]
    peptides = sorted({t["peptide"] for window in windows if window["native_event_review_tier"] == "PASS"
                       for t in target_peptides(window)})
    automaton = ahocorasick.Automaton()
    by_length = defaultdict(set)
    for peptide in peptides:
        automaton.add_word(peptide, peptide)
        by_length[len(peptide)].add(peptide)
    automaton.make_automaton()
    sources_by_peptide = {peptide: set() for peptide in peptides}
    allowed = {f.name for f in dataclasses.fields(SelfReferenceSource)}
    for i, record in enumerate(inventory):
        sequence = record["sequence"]
        hits = {peptide for _, peptide in automaton.iter(sequence)}
        # Independent exhaustive enumeration implements p in sequence literally
        # for every possible query length, including negative matches.
        enumerated = set()
        for length, wanted in by_length.items():
            enumerated.update(sequence[start:start + length] for start in range(len(sequence) - length + 1)
                              if sequence[start:start + length] in wanted)
        assert hits == enumerated
        sources = []
        for source in record["sources"]:
            assert all(value in (None, "", False) for key, value in source.items() if key not in allowed)
            sources.append(SelfReferenceSource(**{k: v for k, v in source.items() if k in allowed}))
        for peptide in hits:
            sources_by_peptide[peptide].update(sources)
        if i % 10000 == 0:
            print("REFERENCE_PROTEINS_SCANNED", i, "/", len(inventory), flush=True)
    result = {}
    for peptide, sources in sources_by_peptide.items():
        ordered = tuple(sorted(sources, key=lambda s: (s.gene_id, s.transcript_id, s.protein_id,
                                                     s.gene_name, s.species, s.transcript_biotype)))
        result[peptide] = dataclasses.asdict(SelfReferenceMatch(peptide=peptide, occurs=bool(ordered),
             antigen_kind="splice", excluded_gene_ids=(), sources=ordered,
             source_provenance_complete=True, genome_release="103", reference_sha256=SNAPSHOT_SHA))
    out = root / "reference-check"
    out.mkdir(exist_ok=True)
    with gzip.open(out / "matches.json.gz", "wt") as handle:
        json.dump(result, handle)
    receipt = dict(status="PASS", inventory_sha256=SNAPSHOT_SHA, reference_protein_fasta_sha256=PROTEIN_SHA,
         reference_gtf_sha256=digest(GTF), reconstruction_sha256=digest(data),
         native_reference_proteins=len(inventory), native_reference_sources=len(associations),
         current_annotation_associations_verified=True, current_protein_sequences_verified=True,
         all_reference_proteins_checked_by_two_exact_matchers=True,
         native_factory_AST_verified=True, native_matching_policy_AST_verified=True,
         peptides=len(peptides), exact_reference_peptides=sum(r["occurs"] for r in result.values()),
         output_sha256=digest(out / "matches.json.gz"), worker_sha256=digest(Path(__file__)))
    (out / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2), flush=True)


if __name__ == "__main__":
    main(Path(sys.argv[1]))
