"""Deterministic, sequence-checked annotations for the known upstream DNA segment.

Uses local, pinned evidence only. Names from reference records are annotations;
they do not change the supplied backbone or establish its circular completeness.
"""
import csv
import hashlib
import json

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation, AfterPosition
from Bio.SeqRecord import SeqRecord


def require(condition, message):
    if not condition:
        raise ValueError(message)


def load_inputs(here):
    sources = json.loads((here / "annotation-sources.json").read_text())
    paths = []
    for key in ("selected_regions", "routing_architecture"):
        entry = sources[key]
        path = here / entry["path"]
        require(hashlib.sha256(path.read_bytes()).hexdigest() == entry["sha256"],
                "Annotation input checksum mismatch: " + entry["path"])
        paths.append(entry["path"])
    with (here / paths[0]).open(newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    selected = {row["region_id"]: row for row in rows}
    require(len(selected) == len(rows) and len(rows) > 0, "Expected unique selected antigen annotation records.")
    architecture = json.loads((here / paths[1]).read_text())
    return sources, selected, architecture, paths


def backbone_location(template, reference):
    dna = reference["reference_dna"]
    require(hashlib.sha256(dna.encode("ascii")).hexdigest() == reference["reference_dna_sha256"],
            "Backbone reference checksum mismatch: " + reference["name"])
    # Explicit recorded differences, rather than permissive similarity thresholds.
    for insertion in sorted(reference.get("insertions", []),
                            key=lambda item: item["offset_0based"], reverse=True):
        at = insertion["offset_0based"]
        require(0 <= at <= len(dna), "Invalid reference insertion coordinate.")
        dna = dna[:at] + insertion["dna"] + dna[at:]
    matches = []
    for strand, query in ((1, dna), (-1, str(Seq(dna).reverse_complement()))):
        matches.extend(SimpleLocation(i, i + len(query), strand=strand)
                       for i in range(len(template) - len(query) + 1)
                       if template.startswith(query, i))
    require(len(matches) == 1, "Expected one complete backbone feature match: " + reference["name"])
    return matches[0]


def _coordinate_record(template, protein, components, cfg, sources, selected, architecture):
    """Annotate the known sequence before the unspecified PolyA boundary."""
    offset = len(cfg["template_context"]["prefix_before_transcribed_sequence_dna"])
    transcript_end = offset + components[-1]["rna_end_exclusive"]
    coding_start = offset + len(cfg["fixed_elements"]["utr5_dna"])
    coding_end = coding_start + len(protein) * 3
    promoter_start = offset - len(cfg["template_context"]["promoter_dna"])
    record = SeqRecord(Seq(template), id="YASH_INITIAL", name="Yash_Initial",
                       description="Yash Design #1 SEC-antigen-MITD construct in supplied RNA Foundry (PL3092) DNA context")
    record.annotations.update(
        molecule_type="DNA", date=sources["record_date"], data_file_division="SYN",
        source="synthetic construct", organism="synthetic construct",
        taxonomy=["artificial sequences"],
        comment=[
            "Local design identifier, not an assigned INSDC accession. " + sources["topology"],
            "Generated from designs/initial/construct/mrna/run.py and genbank.py. "
            "Pinned annotation evidence: annotation-sources.json; selected antigen "
            "records: designs/initial/selection/selected.tsv. Stable region_id and "
            "event_id values are preserved in standard /note qualifiers.",
            "RNA Foundry backbone sequence source: supplied PL3092 example. "
            "Cap chemistry and nucleotide modifications are unspecified. "
            "Feature coordinates describe the designed sequence, not measured "
            "cleavage or antigen presentation.",
        ])

    def add(kind, name, start, end, strand=1, **qualifiers):
        feature = SeqFeature(SimpleLocation(start, end, strand=strand), type=kind,
                             qualifiers={**({"standard_name": [name]} if name else {}), **{
                                 key: value if isinstance(value, list) else [str(value)]
                                 for key, value in qualifiers.items()}})
        record.features.append(feature)
        return feature

    add("source", None, 0, len(template),
        organism="synthetic construct", mol_type="other DNA",
        note="Engineered polyantigen sequence in supplied RNA Foundry/PL3092 context; topology unspecified.")
    backbone = []
    for reference in sources["backbone_features"]:
        location = backbone_location(template, reference)
        a, b = int(location.start), int(location.end)
        require(b <= promoter_start,
                "Backbone annotation overlaps the transcription cassette.")
        notes = [reference["note"],
                 "Reference: " + reference["reference_accession"] + ":" + reference["reference_location"],
                 "Sequence source: " + reference["reference_url"],
                 "Feature naming/boundaries: " + reference["annotation_source_name"],
                 "Full source URLs and reference sequence: annotation-sources.json"]
        extra = {}
        if reference["type"] == "CDS":
            translated = str(location.extract(record.seq).translate(table=11, cds=True))
            require(translated == reference["reference_protein"], "Backbone CDS translation mismatch.")
            extra.update(gene=reference["gene"], product=reference["product"],
                         codon_start="1", transl_table="11", translation=translated)
        add(reference["type"], reference["name"], a, b, location.strand, note=notes, **extra)
        backbone.append((a, b))

    # Explicitly identify remaining context without guessing its biological roles.
    for left, right in ((0, promoter_start),):
        cursor = left
        for a, b in sorted((a, b) for a, b in backbone if left <= a < right):
            require(a >= cursor, "Overlapping backbone feature definitions.")
            if a > cursor:
                add("misc_feature", "Uncharacterized backbone context", cursor, a,
                    note="Preserved supplied DNA; no specific functional annotation assigned.")
            cursor = b
        if cursor < right:
            add("misc_feature", "Uncharacterized backbone context", cursor, right,
                note="Preserved supplied DNA; no specific functional annotation assigned.")

    add("regulatory", "T7 promoter", promoter_start, offset, regulatory_class="promoter",
        note="Promoter sequence and intended transcription start retained from the supplied PL3092 example.")
    add("mRNA", "Intended vaccine transcript", offset, transcript_end,
        product="Yash Design #1 polyantigen mRNA",
        note="Internal known-segment coordinate view; PolyA is unspecified and has no encoded sequence in this view.")
    add("CDS", "SEC-antigens-MITD", coding_start, coding_end + 3,
        product="Yash polyantigen protein with SEC and MITD routing elements",
        translation=protein, codon_start="1", transl_table="1",
        note="Includes the first stop codon. The two additional stop triplets are outside this CDS.")

    seen = set()
    join_number = 0
    for index, part in enumerate(components):
        kind = part["kind"]
        feature_type = "misc_feature"
        notes = ["component_index=" + str(index), "component_kind=" + kind]
        name = part["name"]
        if part["region_id"]:
            row = selected[part["region_id"]]
            require(row["sequence"] == part["aa"] and row["antigen_id"] == part["name"],
                    "Antigen annotation does not match the encoded region.")
            seen.add(part["region_id"])
            category = {"MUTATION": "mutation", "SPLICE": "splice",
                        "CTA/SELF": "CTA/self", "FUSION": "fusion"}[row["source"]]
            gene = row["gene"].replace("--", "::")
            change = (" " + row["protein_change"]) if row["protein_change"] else ""
            name = f'{part["name"]} | {gene}{change} | {category}'
            # A space lets the GenBank writer wrap before a long event ID,
            # instead of inserting whitespace into the identifier itself.
            notes.extend("%s= %s" % (key, row[key]) for key in
                         ("region_id", "event_id", "gene", "source", "protein_change",
                          "context_id", "target_mode", "Caris_flags", "somatic_status") if row.get(key))
            notes.append("Source record: designs/initial/selection/selected.tsv")
            if row["source"] == "SPLICE":
                require(row.get("patient_RNA_sequence_confirmed") == "True"
                        and row.get("upstream_observed_stop_codons_excluded") == "True"
                        and int(row.get("fragment_count", "0")) >= 2,
                        "Splice antigen lacks complete patient-RNA support in a stop-free annotated upstream frame.")
                notes.append("Actual patient RNA bases encode this complete 25-aa region in a GENCODE 37 upstream CDS frame, with no observed intervening stop codon.")
                notes.append("Full-window supporting RNA fragments: " + row["fragment_count"])
                notes.append("Supporting reference transcripts: " + row["supporting_transcripts"])
        elif kind in ("signal_peptide", "mitd"):
            element = "SEC" if kind == "signal_peptide" else "MITD"
            source = architecture["elements"][element]
            require(part["aa"] == source["aa"], "Routing-element sequence differs from its source.")
            name = element + " - Yash's " + source["patient_allele"]
            interval = source["protein_interval_1based_inclusive"]
            notes += ["Patient allele: " + source["patient_allele"],
                      "Native protein residues: %d-%d" % tuple(interval),
                      "Reference database: https://www.ebi.ac.uk/ipd/imgt/hla/",
                      "Reference allele: " + source["reference_allele"],
                      "Reference commit: " + source["reference_commit"],
                      source["sequence_basis"].replace("every encoded base", "every base in the native allele CDS"),
                      "The vaccine coding sequence is codon-optimised for this native protein segment.",
                      "DNA/RNA evidence: designs/initial/construct/allele-evidence/"]
            if kind == "signal_peptide":
                notes.append("26-aa HLA-C segment: annotated 24-aa leader plus the following native CS residues; cleavage has not been measured in this construct.")
        elif kind == "linker":
            join_number += 1
            name = f'Join {join_number:02d}: ' + (part["aa"] or "direct")
            notes.append("Encoded linker: " + part["aa"] if part["aa"] else
                         "Direct join; zero inserted nucleotides. The caret marks a boundary.")
        else:
            feature_type, name, explanation = {
                "5UTR_Kozak": ("5'UTR", "RNA Foundry 5' UTR and upstream Kozak context", "63 nt retained from PL3092; the downstream AUG belongs to the CDS."),
                "stop_codon": ("misc_feature", "Primary stop: TAA", "First in-frame stop; included in the CDS."),
                "additional_stops": ("misc_feature", "Additional stops: TAG TGA", "Two in-frame stop triplets after the CDS."),
                "3UTR": ("3'UTR", "RNA Foundry 3' UTR body", "71-nt body after the three-stop block; followed by unspecified PolyA."),
                "polyA_unspecified": ("misc_feature", "PolyA", "PolyA follows this segment. Length, sequence and format are unspecified. This caret marks an unknown-length boundary, not a zero-length tail."),
            }[kind]
            notes.append(explanation)
        if part["aa"]:
            a = (part["rna_start_0based"] - len(cfg["fixed_elements"]["utr5_dna"])) // 3
            notes.append(f'Encoded protein coordinates: {a + 1}..{a + len(part["aa"])} (1-based inclusive).')
        add(feature_type, name, part["template_start_0based"], part["template_end_exclusive"], note=notes)
    require(seen == set(selected), "Not all selected regions have annotations.")
    record.features.sort(key=lambda f: (f.type != "source", int(f.location.start),
                                        len(f) != 0, -int(f.location.end)))
    return record



def build_record(before, protein, components, cfg, sources, selected, architecture):
    """Return the annotated upstream segment without assigning a 3′ endpoint."""
    record = _coordinate_record(before, protein, components, cfg, sources, selected, architecture)
    offset = len(cfg["template_context"]["prefix_before_transcribed_sequence_dna"])
    record.id = "YASH_INITIAL_5P"
    record.name = "Yash_Initial_5P"
    record.description = "Yash Design #1 known DNA segment before unspecified PolyA; partial construct"
    record.annotations["comment"].append(
        "PARTIAL CONSTRUCT: only the known upstream sequence is exported. "
        "PolyA length and format, any downstream extension and the template "
        "boundary depend on the final architecture and are unspecified. "
        "No complete RNA/template sequence or total length is assigned.")
    for feature in record.features:
        if feature.type == "source":
            feature.qualifiers["note"] = [
                "Known upstream DNA segment only. PolyA and any downstream extension are unspecified."]
        elif feature.type == "mRNA":
            feature.location = SimpleLocation(offset, AfterPosition(len(before)), strand=1)
            feature.qualifiers = {
                "standard_name": ["Known vaccine transcript segment"],
                "product": ["Yash Design #1 partial mRNA"],
                "note": ["Sequence before unspecified PolyA; not a complete transcript."]}
    return [record]


def verify_export(path, records, components, protein, rna_segments, selected):
    """Verify the upstream record and the complete, unchanged coding sequence."""
    exported = list(SeqIO.parse(path, "genbank"))
    require(len(exported) == len(records) == len(rna_segments) == 1,
            "Expected one known upstream DNA record and RNA segment.")
    annotated = {}
    coding_seen = False
    for record_index, (expected, actual) in enumerate(zip(records, exported)):
        require(actual.id == expected.id and actual.seq == expected.seq, "GenBank segment changed.")
        require(len(actual.features) == len(expected.features), "GenBank feature lost on export.")
        for before, feature in zip(expected.features, actual.features):
            require(before.type == feature.type and before.location == feature.location
                    and before.qualifiers == feature.qualifiers, "GenBank annotation changed on round-trip.")
            require(0 <= int(feature.location.start) <= int(feature.location.end) <= len(actual),
                    "GenBank feature lies outside its segment.")
            if feature.type == "CDS":
                translated = str(feature.extract(actual.seq).translate(
                    table=int(feature.qualifiers["transl_table"][0]), cds=True))
                require(translated == feature.qualifiers["translation"][0], "CDS translation mismatch.")
                if feature.qualifiers.get("standard_name") == ["SEC-antigens-MITD"]:
                    require(translated == protein, "Vaccine protein changed.")
                    coding_seen = True
            for note in feature.qualifiers.get("note", []):
                if note.startswith("component_index="):
                    index = int(note.split("=", 1)[1])
                    require(index not in annotated, "Duplicate component annotation.")
                    annotated[index] = (record_index, feature)
            require("region_id" not in feature.qualifiers and "label" not in feature.qualifiers,
                    "Nonstandard GenBank qualifier.")
        source = [f for f in actual.features if f.type == "source"]
        require(len(source) == 1 and len(source[0]) == len(actual), "Missing segment source.")
        transcript = next(f for f in actual.features if f.type == "mRNA")
        require(str(transcript.extract(actual.seq)).replace("T", "U") == rna_segments[record_index],
                "Known RNA/DNA segment mismatch.")
    require(coding_seen, "Missing vaccine CDS.")
    require(set(annotated) == set(range(len(components))), "Missing component annotation.")
    for index, part in enumerate(components):
        segment, feature = annotated[index]
        require(int(feature.location.start) == part["template_start_0based"]
                and int(feature.location.end) == part["template_end_exclusive"],
                "Component segment coordinate mismatch.")
        require(str(feature.extract(exported[segment].seq)) == part["dna"], "Component sequence mismatch.")
        if part["aa"]:
            require(str(feature.extract(exported[segment].seq).translate()) == part["aa"],
                    "Component translation mismatch.")
        if part["region_id"]:
            row = selected[part["region_id"]]
            require(all(key+"= "+row[key] in feature.qualifiers["note"]
                        for key in ("region_id", "gene", "source", "event_id")), "Antigen identity lost.")
    return dict(status="KNOWN_SEGMENTS_VERIFIED", record_count=1,
                feature_count=sum(len(r.features) for r in exported), antigen_regions=len(selected),
                round_trip_annotations_identical=True, component_coordinates_and_sequences_verified=True,
                CDS_translations_verified=True, transcript_segments_match_RNA=True,
                backbone_features_sequence_checked=True, complete_construct=False,
                polyA_status="unspecified", topology="unspecified")
