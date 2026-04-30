"""Somatic-mutation preprocessing.

Builds a (patient x gene) {0,1} matrix of mutation events. Works with
either the long-form Xena `mutect2_snv.tsv.gz` (one row per variant,
columns include `Sample`, `gene`, `effect`/`Variant_Classification`) or
a pre-aggregated patient x gene matrix.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Iterable

import pandas as pd

from ..utils.logging import get_logger
from ._common import keep_primary_tumor, read_tsv, tcga_sample_to_patient

# Effects to keep when `require_nonsynonymous` is true (Xena `effect` values
# follow the MAF Variant_Classification vocabulary).
_NONSYN_EFFECTS = {
    "missense_variant", "missense_mutation", "stop_gained", "nonsense_mutation",
    "frameshift_variant", "frame_shift_del", "frame_shift_ins",
    "inframe_insertion", "inframe_deletion", "in_frame_del", "in_frame_ins",
    "splice_acceptor_variant", "splice_donor_variant", "splice_site",
    "start_lost", "translation_start_site", "stop_lost",
    "protein_altering_variant",
}


def _is_long_format(df: pd.DataFrame) -> bool:
    cols = {c.lower() for c in df.columns}
    return (("sample" in cols) or ("tumor_sample_barcode" in cols)) and (("gene" in cols) or ("hugo_symbol" in cols))


def _from_long(df: pd.DataFrame, top_n: int, require_nonsyn: bool, log) -> pd.DataFrame:
    df = df.rename(columns={c: c.lower() for c in df.columns})
    sample_col = "sample" if "sample" in df.columns else "tumor_sample_barcode"
    gene_col = "gene" if "gene" in df.columns else "hugo_symbol"
    effect_col = None
    for c in ["effect", "variant_classification", "one_consequence", "consequence"]:
        if c in df.columns:
            effect_col = c
            break

    df["patient_id"] = df[sample_col].map(tcga_sample_to_patient)

    if require_nonsyn and effect_col is not None:
        norm = df[effect_col].astype(str).str.lower().str.strip()
        before = len(df)
        df = df[norm.isin(_NONSYN_EFFECTS)]
        log.info(f"Filtered non-synonymous: {before} -> {len(df)} variants")

    df = df.dropna(subset=["patient_id", gene_col]).drop_duplicates(["patient_id", gene_col])
    counts = df.groupby(gene_col).size().sort_values(ascending=False)
    top_genes: Iterable[str] = counts.head(top_n).index
    log.info(f"Top-{top_n} mutated genes selected (max count = {int(counts.iloc[0]) if len(counts) else 0})")

    df = df[df[gene_col].isin(top_genes)]
    mat = (
        pd.crosstab(df["patient_id"], df[gene_col])
        .clip(upper=1)
        .astype("int8")
    )
    mat.index.name = "patient_id"
    return mat


def _from_matrix(df: pd.DataFrame, top_n: int, log) -> pd.DataFrame:
    log.info(f"Treating mutation file as patient x gene matrix ({df.shape})")
    df = df.copy()
    if df.index.name is None:
        df = df.set_index(df.columns[0])
    df.index = [tcga_sample_to_patient(str(s)) for s in df.index]
    df = df.groupby(level=0).max()
    df = (df.fillna(0) > 0).astype("int8")
    if df.shape[1] > top_n:
        keep = df.sum(axis=0).sort_values(ascending=False).head(top_n).index
        df = df.loc[:, keep]
    df.index.name = "patient_id"
    return df


def preprocess_mutations(mutation_path: Path, cfg: Dict[str, Any]) -> pd.DataFrame:
    log = get_logger("preprocess_mutations", log_dir=cfg["paths"]["logs"], level=cfg.get("log_level", "INFO"))
    spec = cfg["preprocessing"]["mutations"]

    df = read_tsv(mutation_path)
    log.info(f"Loaded mutation file: {df.shape}")

    if _is_long_format(df):
        # Filter to tumor samples first
        if "Tumor_Sample_Barcode" in df.columns:
            df = df[df["Tumor_Sample_Barcode"].isin(keep_primary_tumor(df["Tumor_Sample_Barcode"].unique().tolist()))]
        elif "Sample" in df.columns:
            df = df[df["Sample"].isin(keep_primary_tumor(df["Sample"].unique().tolist()))]
        elif "sample" in df.columns:
            df = df[df["sample"].isin(keep_primary_tumor(df["sample"].unique().tolist()))]
        out = _from_long(
            df,
            top_n=int(spec.get("top_n_genes", 300)),
            require_nonsyn=bool(spec.get("require_nonsynonymous", True)),
            log=log,
        )
    else:
        out = _from_matrix(df, top_n=int(spec.get("top_n_genes", 300)), log=log)

    log.info(f"Final mutation matrix: {out.shape}; mean #mutated genes / patient = {out.sum(axis=1).mean():.1f}")
    return out
