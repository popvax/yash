"""Stable local identifiers for the two recorded Design 4 research constructs."""


def variant_for_regions(regions):
    genes = {r["gene"] for r in regions}
    if len(regions) == 16 and "KCTD19" in genes:
        return "MAIN16"
    if len(regions) == 15 and "KCTD19" not in genes:
        return "STRAND15"
    raise ValueError("Expected the recorded MAIN16 or KCTD19-excluded STRAND15 selection")


def export_identifiers(variant):
    if variant not in ("MAIN16", "STRAND15"):
        raise ValueError("Unknown Design 4 construct variant: " + str(variant))
    prefix = "Yash_Design4_" + variant
    return {
        "protein_fasta": prefix + "_RESEARCH_ONLY",
        "cds_fasta": prefix + "_CDS_DNA_including_first_stop",
        "rna_fasta": prefix + "_RNA_before_PolyA_partial_PolyA_unspecified",
        "template_fasta": prefix + "_DNA_before_PolyA_partial_PolyA_unspecified",
        "genbank_locus": "YASH_D4_" + variant,
        "genbank_definition": "Yash Design 4 " + variant + " known DNA segment before unspecified PolyA; partial research construct",
    }


def configure_identifiers(config, variant):
    config.update(name="design-4-" + variant.lower() + "-research-mrna",
                  construct_variant=variant, export_identifiers=export_identifiers(variant))
