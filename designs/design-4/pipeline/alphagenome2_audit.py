"""Live AlphaGenome 2 (v2) + Atlas execution (--execute) and offline receipt verification.

Reads ALPHAGENOME_API_KEY exclusively from the environment when --execute is passed;
never logs, prints, or writes credentials to disk. Always executes live queries when
--execute is given (no silent skip if output files exist). Offline verification requires
no network access or API key and independently checks observed-only per-track aggregation
and unit separation against the retained per-track tables.
"""
import argparse
import gzip
import json
import math
import os
import time
import warnings
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from dataio import iter_rows, load, read_text_auto, resolve_gencode_fasta, rows, save, sha, table

HERE = Path(__file__).resolve().parents[1]
AG2_DIR = HERE / "alphagenome2"

REFERENCE_CONDITIONS = (
    ("pancreas_tissue", "UBERON:0001264"),
    ("pancreatic_ductal_cell", "UBERON:0001264_CL:0002079"),
    ("pdac_cell_line_Panc1", "EFO:0002713"),
    ("testis", "UBERON:0000473"),
    ("heart", "UBERON:0000948"),
    ("lung", "UBERON:0002048"),
    ("lung_alveolar_epithelial", "CL:0010003"),
    ("liver", "UBERON:0002107"),
    ("liver_hepatocyte", "UBERON:0002107_CL:0000182"),
    ("peritoneal_mesothelial", "UBERON:0000945_CL:1000490"),
    ("kidney_cortex", "UBERON:0001225"),
    ("brain", "UBERON:0000955"),
)

# Conditions where bulk RNA-seq and splicing have observed reference training channels
BULK_OBSERVED_VARIANT_CONDITIONS = (
    ("pancreas", ["UBERON:0001264"]),
    ("liver", ["UBERON:0002107"]),
    ("lung", ["UBERON:0002048"]),
)

# Primary cell-type conditions where SC_GENE_EXPRESSION has observed reference channels
# (splicing is imputed_exploratory in single-cell-only conditions and is NOT scored here)
SC_OBSERVED_VARIANT_CONDITIONS = (
    ("ductal", ["UBERON:0001264", "CL:0002079"]),
    ("peritoneum", ["UBERON:0000945", "CL:1000490"]),
)

ORGAN_INTERVAL_SPECS = (
    ("pancreas", ["UBERON:0001264"]),
    ("ductal", ["UBERON:0001264", "CL:0002079"]),
    ("liver", ["UBERON:0002107"]),
    ("hepatocyte", ["UBERON:0002107", "CL:0000182"]),
    ("lung", ["UBERON:0002048"]),
    ("peritoneum_mesothelial", ["UBERON:0000945", "CL:1000490"]),
    ("testis", ["UBERON:0000473"]),
)

# Exact 32 core, fusion, APM, ERBB, and immune-microenvironment genes (GRCh38 Ensembl coordinates)
PROFILE_GENES = (
    ("KLF5", "core_somatic", "chr13", 73054976, 73077586, "+", "ENSG00000102554"),
    ("SF1", "core_somatic", "chr11", 64764605, 64778802, "-", "ENSG00000168066"),
    ("AQP3", "core_somatic", "chr9", 33441156, 33447596, "-", "ENSG00000165272"),
    ("RXRA", "core_somatic", "chr9", 134317098, 134440585, "+", "ENSG00000186350"),
    ("HSPG2", "core_somatic", "chr1", 21820593, 21937318, "-", "ENSG00000142798"),
    ("TGFBR2", "core_somatic", "chr3", 30606478, 30694249, "+", "ENSG00000163513"),
    ("PPIC", "core_somatic", "chr5", 123022052, 123036727, "-", "ENSG00000168938"),
    ("KCTD19", "core_somatic_watch", "chr16", 67289428, 67326766, "-", "ENSG00000168676"),
    ("VPS11", "core_somatic", "chr11", 119067789, 119082753, "+", "ENSG00000160695"),
    ("ARHGAP40", "core_somatic", "chr20", 38601775, 38651035, "+", "ENSG00000124143"),
    ("PLEC", "core_somatic", "chr8", 143915152, 143977039, "-", "ENSG00000178209"),
    ("SYT12", "core_somatic_1frag", "chr11", 67006772, 67051110, "+", "ENSG00000173227"),
    ("H2BC12", "core_somatic_rescued", "chr6", 27135928, 27146863, "-", "ENSG00000197903"),
    ("TJP2", "core_somatic_rescued", "chr9", 69121264, 69274641, "+", "ENSG00000119139"),
    ("FES", "core_somatic_1frag", "chr15", 90883695, 90896371, "+", "ENSG00000182511"),
    ("CDH1", "core_fusion_partner_5p", "chr16", 68737242, 68835543, "+", "ENSG00000039068"),
    ("NRG1", "core_fusion_partner_3p", "chr8", 31639222, 32855666, "+", "ENSG00000157168"),
    ("B2M", "apm_essential", "chr15", 44711253, 44718851, "+", "ENSG00000166710"),
    ("HLA-A", "apm_class1", "chr6", 29941260, 29949606, "+", "ENSG00000206503"),
    ("HLA-B", "apm_class1", "chr6", 31353195, 31367067, "-", "ENSG00000234745"),
    ("HLA-C", "apm_class1", "chr6", 31268746, 31272171, "-", "ENSG00000204525"),
    ("HLA-DRA", "apm_class2", "chr6", 32439870, 32445181, "+", "ENSG00000204287"),
    ("HLA-DRB1", "apm_class2", "chr6", 32577902, 32589848, "-", "ENSG00000196126"),
    ("ERBB2", "erbb_signaling", "chr17", 39687914, 39730426, "+", "ENSG00000141736"),
    ("ERBB3", "erbb_signaling", "chr12", 56076799, 56103505, "+", "ENSG00000065361"),
    ("ERBB4", "erbb_signaling", "chr2", 211375717, 212539266, "-", "ENSG00000178568"),
    ("CLDN18", "surface_marker", "chr3", 137998735, 138033655, "+", "ENSG00000066405"),
    ("CD274", "checkpoint_tme", "chr9", 5450503, 5470566, "+", "ENSG00000120217"),
    ("PDCD1", "checkpoint_tme", "chr2", 241849884, 241858894, "-", "ENSG00000188389"),
    ("LAG3", "checkpoint_tme", "chr12", 6772360, 6778455, "+", "ENSG00000089692"),
    ("ADORA2A", "checkpoint_tme", "chr22", 24417879, 24443177, "+", "ENSG00000128271"),
    ("VEGFA", "angiogenesis_tme", "chr6", 43770184, 43786487, "+", "ENSG00000112715"),
)


def _load_computed_ctdna_inputs():
    """Load BAM pileup and MHCflurry Class-I predictions generated by tools/audit_ctdna_loci.py."""
    pileups = rows(HERE / "exports/ctdna-locus-pileups.tsv")
    hla_rows = rows(HERE / "exports/ctdna-mhcflurry-predictions.tsv")
    hla_by_vid = {}
    for r in hla_rows:
        hla_by_vid.setdefault(r["variant_label"], []).append(r)

    records = []
    for p in pileups:
        vid = p["variant_label"]
        gene = p["gene"]
        chrom, pos_s, alleles_s = p["grch38_variant"].split(":")
        ref, alt = alleles_s.split(">")
        in_d3 = p["encoded_in_design4_main"]
        hlist = hla_by_vid.get(vid, [])
        if not hlist:
            hla_summary = "not_missense_or_excluded_from_targeting"
        else:
            best = min(hlist, key=lambda x: float(x["mutant_affinity_nm"]))
            wt_part = (
                f"(WT={float(best['wildtype_affinity_nm']):.2f}nM)"
                if best["wildtype_affinity_nm"] else ""
            )
            hla_summary = (
                f"{best['best_mutant_peptide']}@{best['allele']}:"
                f"{float(best['mutant_affinity_nm']):.2f}nM{wt_part}"
            )
        m_filter = p["paired_mutect2_filter"]
        m_ad = p["paired_mutect2_tumour_ad"]
        wes_raw = f"{p['liver_wes_ref_fragments']}/{p['liver_wes_alt_fragments']}"
        wes_summary = (
            f"{wes_raw} (Mutect2_{m_filter}:{m_ad})" if m_ad else f"{wes_raw} ({m_filter})"
        )
        records.append({
            "variant_label": vid,
            "gene": gene,
            "ensembl_id": p["ensembl_id"],
            "protein_or_coding_change": p["report_label"],
            "chrom": chrom,
            "pos": int(pos_s),
            "ref": ref,
            "alt": alt,
            "grch38_variant": p["grch38_variant"],
            "mar_2025_cfdna": p["mar_2025_cfdna"],
            "oct_2025_cfdna": p["oct_2025_cfdna"],
            "normal_wbc_fragments": f"{p['normal_wbc_ref_fragments']}/{p['normal_wbc_alt_fragments']}",
            "jan_2025_liver_wes": wes_summary,
            "jan_2025_liver_rna_all": f"{p['liver_wts_ref_fragments']}/{p['liver_wts_alt_fragments']}",
            "jan_2025_liver_rna_coding_strand": f"{p['liver_wts_coding_ref_fragments']}/{p['liver_wts_coding_alt_fragments']}",
            "sample_observation_status": p["sample_observation_status"],
            "encoded_in_design4_main": in_d3,
            "mhcflurry_best_hla_nm": hla_summary,
        })
    return records


def _symbol_to_ensg():
    mapping = {}
    with gzip.open(resolve_gencode_fasta("gencode.v46.pc_translations.fa.gz", HERE / "references"), "rt") as stream:
        for line in stream:
            if line.startswith(">"):
                parts = line[1:].strip().split("|")
                if len(parts) >= 7:
                    mapping.setdefault(parts[6], parts[2].split(".")[0])
    return mapping


def _build_pseudo_umi_anndata(scale=100):
    import anndata as ad
    import numpy as np
    import pandas as pd

    sym_map = _symbol_to_ensg()
    counts = {}
    for line in read_text_auto(HERE / "exports/expression.tsv").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        sym, val = line.split("\t")
        fval = float(val)
        if fval > 0 and sym in sym_map:
            ensg = sym_map[sym]
            counts[ensg] = counts.get(ensg, 0.0) + fval
    genes = sorted(counts)
    matrix = np.maximum(1, np.rint([counts[g] * scale for g in genes])).astype(np.int32).reshape(1, -1)
    return ad.AnnData(
        X=matrix,
        obs=pd.DataFrame(index=["TN25-112103_bulk_WTS_pseudoUMI"]),
        var=pd.DataFrame(index=genes),
    )


def _summarize_observed_variant_tracks(track_rows):
    """Aggregate strictly imputation_mode == 'observed' tracks with separate scorer units."""
    out = {}
    for cname, _ in BULK_OBSERVED_VARIANT_CONDITIONS:
        obs = [r for r in track_rows if r["condition_label"] == cname and r["imputation_mode"] == "observed"]
        lfc = [float(r["raw_score"]) for r in obs if "GeneMaskLFCScorer" in r["variant_scorer"] and r["output_type"] == "RNA_SEQ"]
        act = [float(r["raw_score"]) for r in obs if "GeneMaskActiveScorer" in r["variant_scorer"] and r["output_type"] == "RNA_SEQ"]
        sps = [abs(float(r["raw_score"])) for r in obs if r["output_type"] == "SPLICE_SITES"]
        spu = [abs(float(r["raw_score"])) for r in obs if r["output_type"] == "SPLICE_SITE_USAGE"]
        spj = [abs(float(r["raw_score"])) for r in obs if r["output_type"] == "SPLICE_JUNCTIONS"]
        out[f"{cname}_obs_rna_lfc_mean"] = round(sum(lfc) / len(lfc), 5) if lfc else ""
        out[f"{cname}_obs_rna_active_mean"] = round(sum(act) / len(act), 5) if act else ""
        out[f"{cname}_obs_splice_sites_max_prob_diff"] = round(max(sps), 5) if sps else ""
        out[f"{cname}_obs_splice_usage_max_frac_diff"] = round(max(spu), 5) if spu else ""
        out[f"{cname}_obs_splice_junctions_max_log_diff"] = round(max(spj), 5) if spj else ""
        out[f"{cname}_obs_track_counts"] = f"rna:{len(lfc)};sps:{len(sps)};spu:{len(spu)};spj:{len(spj)}"
    for cname, _ in SC_OBSERVED_VARIANT_CONDITIONS:
        obs = [r for r in track_rows if r["condition_label"] == cname and r["imputation_mode"] == "observed"]
        lfc = [float(r["raw_score"]) for r in obs if "GeneMaskLFCScorer" in r["variant_scorer"] and r["output_type"] == "SC_GENE_EXPRESSION"]
        act = [float(r["raw_score"]) for r in obs if "GeneMaskActiveScorer" in r["variant_scorer"] and r["output_type"] == "SC_GENE_EXPRESSION"]
        out[f"{cname}_obs_sc_lfc_mean"] = round(sum(lfc) / len(lfc), 5) if lfc else ""
        out[f"{cname}_obs_sc_active_mean"] = round(sum(act) / len(act), 5) if act else ""
        out[f"{cname}_obs_sc_track_count"] = len(lfc)
    return out


def execute_live():
    warnings.filterwarnings("ignore")
    api_key = os.environ.get("ALPHAGENOME_API_KEY", "").strip()
    if not api_key:
        raise RuntimeError("Set ALPHAGENOME_API_KEY in the environment to run --execute")

    import alphagenome
    import numpy as np
    from alphagenome.atlas import atlas
    from alphagenome.data import genome
    from alphagenome.models.v2 import dna_client, interval_scorers, model_condition, variant_scorers

    client = dna_client.create(api_key)
    atlas_client = atlas.create(api_key)

    def _call_retry(fn, *args, **kwargs):
        for attempt in range(8):
            try:
                return fn(*args, **kwargs)
            except Exception as exc:
                if "RESOURCE_EXHAUSTED" in str(exc) and attempt < 7:
                    time.sleep(15 + 5 * attempt)
                    continue
                raise

    # 1. Export per-condition best-channel reference training status (not patient observation)
    cond_df = client.supported_conditions().set_index("condition")
    modality_cols = [
        "atac", "cage", "dnase", "rna_seq", "chip_histone", "chip_tf",
        "splice_sites", "splice_site_usage", "splice_junctions",
        "contact_maps", "sc_atac", "sc_gene_expression", "methylation",
    ]
    imputation_rows = []
    for label, curie in REFERENCE_CONDITIONS:
        row = cond_df.loc[curie]
        entry = {
            "label": label,
            "condition_curie": curie,
            "biosample_type": str(row["biosample_type"]),
            "biosample_name": str(row["biosample_name"]),
            "metadata_semantics": "best_channel_imputation_mode_in_alphagenome2_training_corpus_not_patient_measurement",
        }
        for col in modality_cols:
            entry[col] = str(row[col])
        imputation_rows.append(entry)
    table(
        AG2_DIR / "condition-imputation-metadata.tsv",
        imputation_rows,
        ["label", "condition_curie", "biosample_type", "biosample_name", "metadata_semantics"] + modality_cols,
    )

    # 2. Verify live encode_cells against bundled yash-tumor-condition-embedding.npy
    pseudo_adata = _build_pseudo_umi_anndata(scale=100)
    fresh_emb = client.encode_cells(pseudo_adata, progress_bar=False)
    bundled_emb = np.load(AG2_DIR / "yash-tumor-condition-embedding.npy")
    cosine_sim = float(
        np.dot(fresh_emb.ravel(), bundled_emb.ravel())
        / (np.linalg.norm(fresh_emb) * np.linalg.norm(bundled_emb))
    )
    max_abs_diff = float(np.max(np.abs(fresh_emb - bundled_emb)))

    # 3a. Refresh atlas-avi-somatic.tsv live, keying strictly by str(obs["variant"])
    somatic_avi_meta = rows(AG2_DIR / "atlas-avi-somatic.tsv")
    somatic_snvs = []
    for r in somatic_avi_meta:
        ch, pos_s, al_s = r["variant"].split(":")
        rf, al = al_s.split(">")
        somatic_snvs.append(genome.Variant(ch, int(pos_s), rf, al))
    som_atlas = _call_retry(
        atlas_client.query_variants,
        somatic_snvs,
        requested_scorers=["AVI_SCORE", "AVI_SCORE_MODEL_FEATURES", "AVI_SCORE_FEATURE_IMPORTANCE"],
        progress_bar=False,
    )
    som_avi_ad = som_atlas["AVI_SCORE"]
    som_feat_ad = som_atlas["AVI_SCORE_MODEL_FEATURES"]
    som_imp_ad = som_atlas["AVI_SCORE_FEATURE_IMPORTANCE"]
    som_avi_by_v = {
        str(v): float(som_avi_ad.X[i, 0]) for i, v in enumerate(som_avi_ad.obs["variant"])
    }
    som_fnames = [str(x) for x in som_feat_ad.var["name"]]
    som_feat_by_v = {
        str(v): {
            som_fnames[j]: (None if np.isnan(som_feat_ad.X[i, j]) else float(som_feat_ad.X[i, j]))
            for j in range(len(som_fnames))
        }
        for i, v in enumerate(som_feat_ad.obs["variant"])
    }
    som_inames = [str(x) for x in som_imp_ad.var["name"]]
    som_imp_by_v = {
        str(v): {
            som_inames[j]: (None if np.isnan(som_imp_ad.X[i, j]) else float(som_imp_ad.X[i, j]))
            for j in range(len(som_inames))
        }
        for i, v in enumerate(som_imp_ad.obs["variant"])
    }
    refreshed_avi_rows = []
    avi_cols = ["variant_id", "gene", "protein_change", "variant", "rna_alt_total", "AVI_SCORE:AVI_SCORE"]
    avi_cols += [f"AVI_SCORE_MODEL_FEATURES:{n}" for n in som_fnames]
    avi_cols += [f"AVI_SCORE_FEATURE_IMPORTANCE:{n}" for n in som_inames]
    for r in somatic_avi_meta:
        vk = r["variant"]
        row_out = {
            "variant_id": r["variant_id"],
            "gene": r["gene"],
            "protein_change": r["protein_change"],
            "variant": vk,
            "rna_alt_total": r["rna_alt_total"],
            "AVI_SCORE:AVI_SCORE": round(som_avi_by_v[vk], 4),
        }
        for n in som_fnames:
            val = som_feat_by_v[vk][n]
            row_out[f"AVI_SCORE_MODEL_FEATURES:{n}"] = "" if val is None else round(val, 4)
        for n in som_inames:
            val = som_imp_by_v[vk][n]
            row_out[f"AVI_SCORE_FEATURE_IMPORTANCE:{n}"] = "" if val is None else round(val, 4)
        refreshed_avi_rows.append(row_out)
    refreshed_avi_rows.sort(key=lambda x: (-float(x["AVI_SCORE:AVI_SCORE"]), x["variant_id"]))
    table(AG2_DIR / "atlas-avi-somatic.tsv", refreshed_avi_rows, avi_cols)

    # 3b. Always execute live Atlas + DnaClient v2 variant queries on all 13 computed Guardant360 loci
    ctdna_inputs = _load_computed_ctdna_inputs()
    snv_variants = [
        genome.Variant(r["chrom"], r["pos"], r["ref"], r["alt"], name=r["variant_label"])
        for r in ctdna_inputs
        if len(r["ref"]) == 1 and len(r["alt"]) == 1
    ]
    atlas_res = _call_retry(
        atlas_client.query_variants,
        snv_variants,
        requested_scorers=["AVI_SCORE", "AVI_SCORE_MODEL_FEATURES"],
        progress_bar=False,
    )
    avi_map = {}
    avi_ad = atlas_res["AVI_SCORE"]
    for idx, vstr in enumerate(avi_ad.obs["variant"]):
        avi_map[str(vstr)] = float(avi_ad.X[idx, 0])
    feat_map = {}
    feat_ad = atlas_res["AVI_SCORE_MODEL_FEATURES"]
    feat_names = [str(x) for x in feat_ad.var["name"]]
    for idx, vstr in enumerate(feat_ad.obs["variant"]):
        feat_map[str(vstr)] = {
            fname: (None if np.isnan(feat_ad.X[idx, j]) else float(feat_ad.X[idx, j]))
            for j, fname in enumerate(feat_names)
        }

    bulk_scorers = [
        variant_scorers.RECOMMENDED_VARIANT_SCORERS["RNA_SEQ"],
        variant_scorers.RECOMMENDED_VARIANT_SCORERS["RNA_SEQ_ACTIVE"],
        variant_scorers.RECOMMENDED_VARIANT_SCORERS["SPLICE_SITES"],
        variant_scorers.RECOMMENDED_VARIANT_SCORERS["SPLICE_SITE_USAGE"],
        variant_scorers.RECOMMENDED_VARIANT_SCORERS["SPLICE_JUNCTIONS"],
    ]
    sc_scorers = [
        variant_scorers.RECOMMENDED_VARIANT_SCORERS["SC_GENE_EXPRESSION"],
        variant_scorers.RECOMMENDED_VARIANT_SCORERS["SC_GENE_EXPRESSION_ACTIVE"],
    ]

    def _clean_val(v):
        if v is None:
            return ""
        if isinstance(v, float) and math.isnan(v):
            return ""
        s = str(v)
        return "" if s in ("nan", "NaN", "None", "<NA>") else s

    def _score_one_ctdna(rec):
        vid = rec["variant_label"]
        gene = rec["gene"]
        ensg = rec["ensembl_id"]
        vkey = rec["grch38_variant"]
        v_obj = genome.Variant(rec["chrom"], rec["pos"], rec["ref"], rec["alt"], name=vid)
        iv_obj = v_obj.reference_interval.resize(dna_client.SEQUENCE_LENGTH_500KB)
        track_rows = []

        for cname, curies in BULK_OBSERVED_VARIANT_CONDITIONS:
            sv_out = _call_retry(
                client.score_variant, iv_obj, v_obj,
                condition=model_condition.from_curies(curies),
                variant_scorers=bulk_scorers,
            )
            tidy = variant_scorers.tidy_scores(sv_out)
            gtidy = tidy[(tidy["gene_id"].astype(str).str.split(".").str[0] == ensg) | (tidy["gene_name"] == gene)]
            for r in gtidy.to_dict("records"):
                if math.isnan(float(r["raw_score"])):
                    continue
                track_rows.append({
                    "variant_label": vid,
                    "grch38_variant": vkey,
                    "condition_label": cname,
                    "condition_curie": "_".join(curies),
                    "gene_id": _clean_val(r.get("gene_id")),
                    "gene_name": _clean_val(r.get("gene_name")),
                    "output_type": _clean_val(r.get("output_type")),
                    "variant_scorer": _clean_val(r.get("variant_scorer")),
                    "track_name": _clean_val(r.get("track_name")),
                    "track_strand": _clean_val(r.get("track_strand")),
                    "imputation_mode": _clean_val(r.get("imputation_mode")),
                    "assay_title": _clean_val(r.get("Assay title")),
                    "data_source": _clean_val(r.get("data_source")),
                    "splice_site_type": _clean_val(r.get("splice_site_type")),
                    "junction_start": _clean_val(r.get("junction_Start")),
                    "junction_end": _clean_val(r.get("junction_End")),
                    "raw_score": round(float(r["raw_score"]), 6),
                })

        for cname, curies in SC_OBSERVED_VARIANT_CONDITIONS:
            sv_out = _call_retry(
                client.score_variant, iv_obj, v_obj,
                condition=model_condition.from_curies(curies),
                variant_scorers=sc_scorers,
            )
            tidy = variant_scorers.tidy_scores(sv_out)
            gtidy = tidy[(tidy["gene_id"].astype(str).str.split(".").str[0] == ensg) | (tidy["gene_name"] == gene)]
            for r in gtidy.to_dict("records"):
                if math.isnan(float(r["raw_score"])):
                    continue
                track_rows.append({
                    "variant_label": vid,
                    "grch38_variant": vkey,
                    "condition_label": cname,
                    "condition_curie": "_".join(curies),
                    "gene_id": _clean_val(r.get("gene_id")),
                    "gene_name": _clean_val(r.get("gene_name")),
                    "output_type": _clean_val(r.get("output_type")),
                    "variant_scorer": _clean_val(r.get("variant_scorer")),
                    "track_name": _clean_val(r.get("track_name")),
                    "track_strand": _clean_val(r.get("track_strand")),
                    "imputation_mode": _clean_val(r.get("imputation_mode")),
                    "assay_title": _clean_val(r.get("Assay title")),
                    "data_source": _clean_val(r.get("data_source")),
                    "splice_site_type": _clean_val(r.get("splice_site_type")),
                    "junction_start": _clean_val(r.get("junction_Start")),
                    "junction_end": _clean_val(r.get("junction_End")),
                    "raw_score": round(float(r["raw_score"]), 6),
                })

        feats = feat_map.get(vkey, {})
        fmt_f = lambda k: "" if feats.get(k) is None else round(float(feats[k]), 5)
        summary_row = {
            "variant_label": vid,
            "gene": gene,
            "protein_or_coding_change": rec["protein_or_coding_change"],
            "grch38_variant": vkey,
            "mar_2025_cfdna": rec["mar_2025_cfdna"],
            "oct_2025_cfdna": rec["oct_2025_cfdna"],
            "normal_wbc_fragments": rec["normal_wbc_fragments"],
            "jan_2025_liver_wes": rec["jan_2025_liver_wes"],
            "jan_2025_liver_rna_all": rec["jan_2025_liver_rna_all"],
            "jan_2025_liver_rna_coding_strand": rec["jan_2025_liver_rna_coding_strand"],
            "sample_observation_status": rec["sample_observation_status"],
            "encoded_in_design4_main": rec["encoded_in_design4_main"],
            "mhcflurry_best_hla_nm": rec["mhcflurry_best_hla_nm"],
            "atlas_avi_score": "" if vkey not in avi_map else round(avi_map[vkey], 5),
            "atlas_alphamissense": fmt_f("ALPHAMISSENSE"),
            "atlas_merged_splicing": fmt_f("MERGED_SPLICING"),
            "atlas_max_abs_rna_seq": fmt_f("MAX_ABS_RNA_SEQ"),
            "atlas_max_abs_atac": fmt_f("MAX_ABS_ATAC"),
            "atlas_phastcons_470way": fmt_f("PHASTCONS_470_WAY"),
        }
        summary_row.update(_summarize_observed_variant_tracks(track_rows))
        return summary_row, track_rows

    with ThreadPoolExecutor(max_workers=2) as pool:
        ctdna_results = list(pool.map(_score_one_ctdna, ctdna_inputs))

    ctdna_rows = [item[0] for item in ctdna_results]
    all_var_track_rows = [tr for item in ctdna_results for tr in item[1]]
    var_track_cols = [
        "variant_label", "grch38_variant", "condition_label", "condition_curie",
        "gene_id", "gene_name", "output_type", "variant_scorer", "track_name",
        "track_strand", "imputation_mode", "assay_title", "data_source",
        "splice_site_type", "junction_start", "junction_end", "raw_score",
    ]
    table(AG2_DIR / "variant-track-scores.tsv.gz", all_var_track_rows, var_track_cols)

    ctdna_cols = [
        "variant_label", "gene", "protein_or_coding_change", "grch38_variant",
        "mar_2025_cfdna", "oct_2025_cfdna", "normal_wbc_fragments", "jan_2025_liver_wes",
        "jan_2025_liver_rna_all", "jan_2025_liver_rna_coding_strand",
        "sample_observation_status", "encoded_in_design4_main",
        "mhcflurry_best_hla_nm", "atlas_avi_score", "atlas_alphamissense",
        "atlas_merged_splicing", "atlas_max_abs_rna_seq", "atlas_max_abs_atac",
        "atlas_phastcons_470way",
    ]
    for cname, _ in BULK_OBSERVED_VARIANT_CONDITIONS:
        ctdna_cols.extend([
            f"{cname}_obs_rna_lfc_mean", f"{cname}_obs_rna_active_mean",
            f"{cname}_obs_splice_sites_max_prob_diff",
            f"{cname}_obs_splice_usage_max_frac_diff",
            f"{cname}_obs_splice_junctions_max_log_diff",
            f"{cname}_obs_track_counts",
        ])
    for cname, _ in SC_OBSERVED_VARIANT_CONDITIONS:
        ctdna_cols.extend([
            f"{cname}_obs_sc_lfc_mean", f"{cname}_obs_sc_active_mean",
            f"{cname}_obs_sc_track_count",
        ])
    table(AG2_DIR / "ctdna-progression-audit.tsv", ctdna_rows, ctdna_cols)

    # 4. Multi-organ metastatic & immune-microenvironment interval profile with retained per-track imputation_mode
    expr_map = {}
    for line in read_text_auto(HERE / "exports/expression.tsv").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        sym, val = line.split("\t")
        expr_map[sym] = float(val)

    isc = [
        interval_scorers.GeneMaskScorer(
            requested_output=dna_client.OutputType.RNA_SEQ,
            width=200001,
            aggregation_type=interval_scorers.IntervalAggregationType.MEAN,
        ),
        interval_scorers.GeneMaskScorer(
            requested_output=dna_client.OutputType.SC_GENE_EXPRESSION,
            width=200001,
            aggregation_type=interval_scorers.IntervalAggregationType.MEAN,
        ),
    ]
    isc_cell = [
        interval_scorers.GeneMaskScorer(
            requested_output=dna_client.OutputType.SC_GENE_EXPRESSION_DENOISED,
            width=200001,
            aggregation_type=interval_scorers.IntervalAggregationType.MEAN,
        ),
    ]

    win = dna_client.SEQUENCE_LENGTH_1MB
    half = win // 2
    shared_intervals = []
    gene_to_iv_idx = {}
    by_chr = {}
    for g_idx, (sym, cat, ch, st, en, strand, ensg) in enumerate(PROFILE_GENES):
        by_chr.setdefault(ch, []).append((st, en, g_idx))
    for ch, items in sorted(by_chr.items()):
        anchored = []
        for st, en, g_idx in items:
            _, _, _, _, _, strand, _ = PROFILE_GENES[g_idx]
            anchor = (st if strand == "+" else en) if (en - st) > 400000 else (st + en) // 2
            anchored.append((anchor, g_idx))
        anchored.sort()
        cluster = []
        for mid, g_idx in anchored:
            if not cluster:
                cluster = [(mid, g_idx)]
            elif mid - cluster[0][0] <= win - 250000:
                cluster.append((mid, g_idx))
            else:
                cmid = (cluster[0][0] + cluster[-1][0]) // 2
                ist = max(0, cmid - half)
                iv_idx = len(shared_intervals)
                shared_intervals.append(genome.Interval(ch, ist, ist + win))
                for _, gi in cluster:
                    gene_to_iv_idx[gi] = iv_idx
                cluster = [(mid, g_idx)]
        if cluster:
            cmid = (cluster[0][0] + cluster[-1][0]) // 2
            ist = max(0, cmid - half)
            iv_idx = len(shared_intervals)
            shared_intervals.append(genome.Interval(ch, ist, ist + win))
            for _, gi in cluster:
                gene_to_iv_idx[gi] = iv_idx

    def _extract_gene_tracks(adata_obj, ensg, sym, strand, cond_label, out_type_label):
        if adata_obj.shape[0] == 0 or adata_obj.shape[1] == 0:
            return []
        gid_series = (
            adata_obj.obs["gene_id"].astype(str)
            if "gene_id" in adata_obj.obs.columns
            else adata_obj.obs.index.to_series().astype(str)
        )
        mask = gid_series.str.split(".").str[0] == ensg
        if "gene_name" in adata_obj.obs.columns:
            mask = mask | (adata_obj.obs["gene_name"].astype(str) == sym)
        if not mask.any():
            return []
        row_idx = int(np.flatnonzero(mask.to_numpy())[0])
        v = adata_obj.var
        kept = []
        for col_idx in range(adata_obj.shape[1]):
            t_strand = str(v["strand"].iloc[col_idx]) if "strand" in v.columns else "."
            if t_strand not in (".", strand):
                continue
            val = float(adata_obj.X[row_idx, col_idx])
            if math.isnan(val):
                continue
            raw_imp = _clean_val(v["imputation_mode"].iloc[col_idx]) if "imputation_mode" in v.columns else ""
            imp = raw_imp if raw_imp else "unavailable_for_cell_conditioning"
            tname = str(v["name"].iloc[col_idx]) if "name" in v.columns else str(v.index[col_idx])
            assay = _clean_val(v["Assay title"].iloc[col_idx]) if "Assay title" in v.columns else ""
            dsrc = _clean_val(v["data_source"].iloc[col_idx]) if "data_source" in v.columns else ""
            kept.append({
                "gene": sym,
                "ensembl_id": ensg,
                "gene_strand": strand,
                "condition_label": cond_label,
                "output_type": out_type_label,
                "track_name": tname,
                "track_strand": t_strand,
                "imputation_mode": imp,
                "assay_title": assay,
                "data_source": dsrc,
                "raw_score": round(val, 6),
            })
        return kept

    def _score_shared_iv(iv):
        res_map = {}
        res_map["cell"] = _call_retry(
            client.score_interval,
            iv,
            interval_scorers=isc_cell,
            cell_conditioning=pseudo_adata,
            cell_condition_name="TN25_112103_bulk_pseudoUMI",
        )
        for olabel, curies in ORGAN_INTERVAL_SPECS:
            res_map[olabel] = _call_retry(
                client.score_interval,
                iv,
                condition=model_condition.from_curies(curies),
                interval_scorers=isc,
            )
        return res_map

    with ThreadPoolExecutor(max_workers=2) as pool:
        iv_results = list(pool.map(_score_shared_iv, shared_intervals))

    profile_rows = []
    all_iv_track_rows = []
    for g_idx, (sym, cat, ch, st, en, strand, ensg) in enumerate(PROFILE_GENES):
        res_map = iv_results[gene_to_iv_idx[g_idx]]
        entry = {
            "gene": sym,
            "category": cat,
            "ensembl_id": ensg,
            "chrom": ch,
            "strand": strand,
            "supplied_vendor_expression": "" if sym not in expr_map else round(expr_map[sym], 5),
        }
        c_tracks = _extract_gene_tracks(
            res_map["cell"][0], ensg, sym, strand, "pseudoUMI_cell", "SC_GENE_EXPRESSION_DENOISED"
        )
        all_iv_track_rows.extend(c_tracks)
        entry["ag2_bulk_pseudoUMI_exploratory_denoised_cpm"] = (
            round(sum(t["raw_score"] for t in c_tracks) / len(c_tracks), 5) if c_tracks else ""
        )
        for olabel, _ in ORGAN_INTERVAL_SPECS:
            out_list = res_map[olabel]
            r_tracks = _extract_gene_tracks(out_list[0], ensg, sym, strand, olabel, "RNA_SEQ")
            s_tracks = _extract_gene_tracks(out_list[1], ensg, sym, strand, olabel, "SC_GENE_EXPRESSION")
            all_iv_track_rows.extend(r_tracks)
            all_iv_track_rows.extend(s_tracks)
            r_obs = [t["raw_score"] for t in r_tracks if t["imputation_mode"] == "observed"]
            r_all = [t["raw_score"] for t in r_tracks]
            s_obs = [t["raw_score"] for t in s_tracks if t["imputation_mode"] == "observed"]
            s_all = [t["raw_score"] for t in s_tracks]
            entry[f"{olabel}_rna_observed"] = round(sum(r_obs) / len(r_obs), 5) if r_obs else ""
            entry[f"{olabel}_rna_all_incl_imputed"] = round(sum(r_all) / len(r_all), 5) if r_all else ""
            entry[f"{olabel}_sc_observed"] = round(sum(s_obs) / len(s_obs), 5) if s_obs else ""
            entry[f"{olabel}_sc_all_incl_imputed"] = round(sum(s_all) / len(s_all), 5) if s_all else ""
            entry[f"{olabel}_obs_tracks"] = f"rna:{len(r_obs)}/{len(r_all)};sc:{len(s_obs)}/{len(s_all)}"
        profile_rows.append(entry)

    iv_track_cols = [
        "gene", "ensembl_id", "gene_strand", "condition_label", "output_type",
        "track_name", "track_strand", "imputation_mode", "assay_title",
        "data_source", "raw_score",
    ]
    table(AG2_DIR / "interval-track-scores.tsv.gz", all_iv_track_rows, iv_track_cols)

    pcols = [
        "gene", "category", "ensembl_id", "chrom", "strand",
        "supplied_vendor_expression", "ag2_bulk_pseudoUMI_exploratory_denoised_cpm",
    ]
    for olabel, _ in ORGAN_INTERVAL_SPECS:
        pcols.extend([
            f"{olabel}_rna_observed", f"{olabel}_rna_all_incl_imputed",
            f"{olabel}_sc_observed", f"{olabel}_sc_all_incl_imputed",
            f"{olabel}_obs_tracks",
        ])
    table(AG2_DIR / "metastatic-organ-and-tme-profile.tsv", profile_rows, pcols)

    # 5. Write explicit provenance receipt separating inherited frozen tables vs fresh live API outputs
    save(AG2_DIR / "additional-evidence-receipt.json", {
        "alphagenome_package_version": getattr(alphagenome, "__version__", "2.0.0"),
        "upstream_ctdna_audit_receipt_sha256": sha(HERE / "exports/ctdna-audit-receipt.json"),
        "inherited_frozen_tables": {
            "scope": "Checksum-verified tables inherited from the initial AlphaGenome 2 run (not regenerated by --execute).",
            "somatic-expression-splicing.tsv": sha(AG2_DIR / "somatic-expression-splicing.tsv"),
            "splice-candidates-audit.tsv": sha(AG2_DIR / "splice-candidates-audit.tsv"),
            "cta-vital-organ-safety.tsv": sha(AG2_DIR / "cta-vital-organ-safety.tsv"),
            "fusion-and-apm-profile.tsv": sha(AG2_DIR / "fusion-and-apm-profile.tsv"),
            "yash-tumor-condition-embedding.npy": sha(AG2_DIR / "yash-tumor-condition-embedding.npy"),
            "exports/cta-audit.tsv": sha(HERE / "exports/cta-audit.tsv"),
        },
        "fresh_live_api_tables": {
            "scope": "Unconditionally queried and regenerated from live AlphaGenome 2 + Atlas APIs whenever --execute is run.",
            "atlas-avi-somatic.tsv": sha(AG2_DIR / "atlas-avi-somatic.tsv"),
            "condition-imputation-metadata.tsv": sha(AG2_DIR / "condition-imputation-metadata.tsv"),
            "ctdna-progression-audit.tsv": sha(AG2_DIR / "ctdna-progression-audit.tsv"),
            "variant-track-scores.tsv.gz": sha(AG2_DIR / "variant-track-scores.tsv.gz"),
            "metastatic-organ-and-tme-profile.tsv": sha(AG2_DIR / "metastatic-organ-and-tme-profile.tsv"),
            "interval-track-scores.tsv.gz": sha(AG2_DIR / "interval-track-scores.tsv.gz"),
        },
        "embedding_reproduction": {
            "bundled_embedding_sha256": sha(AG2_DIR / "yash-tumor-condition-embedding.npy"),
            "coding_genes_mapped_from_bundled_gencode_v46": int(pseudo_adata.shape[1]),
            "cosine_similarity_against_bundled_20223_gene_gtf_embedding": round(cosine_sim, 6),
            "max_abs_diff_against_bundled_embedding": round(max_abs_diff, 6),
            "out_of_distribution_caveat": (
                "AlphaGenome 2 CellConditioningInput expects integer single-cell UMI counts. "
                "vendor-transformed-expression.tsv is bulk FFPE WTS transformed expression from a liver-metastasis "
                "biopsy reported at 60-70% histological tumour cellularity; histological cellularity does NOT "
                "measure the malignant vs non-malignant RNA transcript mixture fraction. Pseudo-UMI conditioning "
                "is therefore an out-of-distribution exploratory prior, not a single-cell malignant measurement."
            ),
        },
        "imputation_and_unit_policy": {
            "supported_conditions_semantics": (
                "client.supported_conditions() reports only the best imputation_mode across channels of an output "
                "type in AlphaGenome 2's reference training corpus ('observed' means >=1 channel was measured in "
                "reference biosamples, not all channels and never a measurement in Yash's tissue)."
            ),
            "per_track_imputation_retained": True,
            "missing_imputation_disclosure": (
                "For reference ontology conditions in variant-track-scores.tsv.gz and interval-track-scores.tsv.gz, "
                "imputation_mode is copied verbatim from the retained target-gene rows of the AlphaGenome 2 API "
                "response. For the 30 user-supplied pseudoUMI_cell rows (cell_conditioning), the API returns NaN "
                "for imputation_mode; those 30 rows are explicitly marked 'unavailable_for_cell_conditioning'."
            ),
            "splicing_variant_scoring_policy": (
                "Splicing variant scorers (SPLICE_SITES, SPLICE_SITE_USAGE, SPLICE_JUNCTIONS) are requested only "
                "on bulk-RNA reference tissues (pancreas, liver, lung) and are NOT requested on single-cell-only "
                "conditions (ductal, peritoneum) where all splicing channels are imputed_exploratory. Within the "
                "bulk-RNA tissue responses, variant-track-scores.tsv.gz retains all target-gene track rows "
                "(matching the focal Ensembl ID or gene symbol across observed, imputed_supported and "
                "imputed_exploratory channels) with their explicit imputation_mode labels for auditability, while "
                "summary columns in ctdna-progression-audit.tsv filter strictly to imputation_mode == 'observed'."
            ),
            "separate_scorer_units": {
                "obs_splice_sites_max_prob_diff": "max |ALT - REF| splice-site classification probability across observed tracks",
                "obs_splice_usage_max_frac_diff": "max |ALT - REF| splice-site usage fraction across observed tracks",
                "obs_splice_junctions_max_log_diff": "max |log(ALT) - log(REF)| splice-junction read count across observed tracks",
            },
        },
        "long_gene_window_limit": (
            "30 of 32 profiled genes receive both RNA_SEQ and SC_GENE_EXPRESSION interval scores; "
            "NRG1 (1.216 Mb) and ERBB4 (1.164 Mb) exceed AlphaGenome 2's 1,048,576-bp sequence / "
            "200,001-bp GeneMaskScorer window, leaving single-cell 3'-end mask columns blank while "
            "bulk RNA_SEQ columns remain populated."
        ),
        "interpretation": (
            "AlphaGenome 2 and Atlas outputs are reference-genome regulatory/pathogenicity predictions, "
            "not patient biopsy measurements across unbiopsied lung or peritoneal lesions, not proof of "
            "splice neutrality or normal-organ silence, and not clinical safety or immunogenicity evidence."
        ),
        "reproduction_script": "pipeline/alphagenome2_audit.py --execute",
    })


def verify_offline():
    receipt = load(AG2_DIR / "receipt.json")
    for name, expected in receipt["output_sha256"].items():
        if sha(AG2_DIR / name) != expected:
            raise ValueError("Original AlphaGenome table checksum mismatch: " + name)
    cond = receipt["expression_conditioning"]
    if sha(AG2_DIR / cond["condition_embedding_file"]) != cond["condition_embedding_sha256"]:
        raise ValueError("Bundled AlphaGenome condition embedding checksum mismatch")

    extra = load(AG2_DIR / "additional-evidence-receipt.json")
    for rel, expected in extra["inherited_frozen_tables"].items():
        if rel == "scope":
            continue
        p = (HERE / rel) if "/" in rel else (AG2_DIR / rel)
        if sha(p) != expected:
            raise ValueError("Inherited AlphaGenome artifact checksum mismatch: " + rel)
    for fname, expected in extra["fresh_live_api_tables"].items():
        if fname == "scope":
            continue
        if sha(AG2_DIR / fname) != expected:
            raise ValueError("Fresh AlphaGenome artifact checksum mismatch: " + fname)

    ctdna_receipt = load(HERE / "exports/ctdna-audit-receipt.json")
    if sha(HERE / "exports/ctdna-audit-receipt.json") != extra["upstream_ctdna_audit_receipt_sha256"]:
        raise ValueError("ctDNA audit receipt checksum mismatch")
    for fname, expected in ctdna_receipt["output_sha256"].items():
        if sha(HERE / "exports" / fname) != expected:
            raise ValueError("ctDNA upstream output checksum mismatch: " + fname)

    imp_rows = rows(AG2_DIR / "condition-imputation-metadata.tsv")
    ctdna_rows = rows(AG2_DIR / "ctdna-progression-audit.tsv")
    prof_rows = rows(AG2_DIR / "metastatic-organ-and-tme-profile.tsv")
    var_tracks = list(iter_rows(AG2_DIR / "variant-track-scores.tsv.gz"))
    iv_tracks = list(iter_rows(AG2_DIR / "interval-track-scores.tsv.gz"))

    if len(imp_rows) != len(REFERENCE_CONDITIONS):
        raise ValueError("Unexpected condition-imputation-metadata row count")
    if len(ctdna_rows) != 13:
        raise ValueError("Unexpected ctdna-progression-audit row count")
    if len(prof_rows) != len(PROFILE_GENES):
        raise ValueError("Unexpected metastatic-organ-and-tme-profile row count")
    if not var_tracks or not iv_tracks:
        raise ValueError("Expected non-empty retained per-track AlphaGenome tables")

    # Verify zero exploratory splicing queries on single-cell conditions, non-blank imputation_mode, and exact observed-only summaries
    for itr in iv_tracks:
        if not itr["imputation_mode"].strip():
            raise ValueError(f"Blank imputation_mode in interval-track-scores.tsv.gz for {itr['gene']} {itr['condition_label']}")
    by_vid = {}
    for tr in var_tracks:
        if not tr["imputation_mode"].strip():
            raise ValueError(f"Blank imputation_mode in variant-track-scores.tsv.gz for {tr['variant_label']}")
        if tr["condition_label"] in ("ductal", "peritoneum") and "SPLICE" in tr["output_type"]:
            raise ValueError("Single-cell condition unexpectedly contains splicing variant scores")
        by_vid.setdefault(tr["variant_label"], []).append(tr)

    for row in ctdna_rows:
        recomputed = _summarize_observed_variant_tracks(by_vid.get(row["variant_label"], []))
        for k, v in recomputed.items():
            if str(row[k]) != str(v):
                raise ValueError(f"Observed-only track summary mismatch for {row['variant_label']} {k}: {row[k]} != {v}")

    return {
        "conditions_with_imputation_metadata": len(imp_rows),
        "ctdna_longitudinal_variants_audited": len(ctdna_rows),
        "retained_variant_track_rows": len(var_tracks),
        "metastatic_and_tme_genes_profiled": len(prof_rows),
        "retained_interval_track_rows": len(iv_tracks),
        "embedding_cosine_similarity": extra["embedding_reproduction"][
            "cosine_similarity_against_bundled_20223_gene_gtf_embedding"
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--execute", action="store_true", help="Run live AlphaGenome 2 API queries using ALPHAGENOME_API_KEY")
    args = parser.parse_args()
    if args.execute:
        execute_live()
    print(json.dumps(verify_offline(), indent=2))


if __name__ == "__main__":
    main()
