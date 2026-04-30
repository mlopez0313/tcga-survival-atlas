"""Download TCGA modalities via TCGAbiolinks in R.

This replaces the previous placeholder GDC path with an R-based downloader
that uses TCGAbiolinks inside a micromamba/conda environment.

The downloader writes simple TSV/TSV.GZ files into `data/raw/` with the same
filenames expected elsewhere in the pipeline whenever practical.
"""
from __future__ import annotations

import os
import subprocess
import textwrap
from pathlib import Path
from typing import Any, Dict

from ..utils.logging import get_logger


def _run_r(script: str, env: Dict[str, str], log) -> None:
    rscript_bin = env.get("GDC_RSCRIPT", "Rscript")
    proc = subprocess.run([rscript_bin, "-"], input=script, capture_output=True, text=True, env=env)
    if proc.stdout:
        log.info(proc.stdout.rstrip())
    if proc.stderr:
        log.warning(proc.stderr.rstrip())
    if proc.returncode != 0:
        raise RuntimeError(f"Rscript failed with exit code {proc.returncode}")


def download_gdc_assets(cfg: Dict[str, Any]) -> Dict[str, Path]:
    log = get_logger("download_gdc", log_dir=cfg["paths"]["logs"], level=cfg.get("log_level", "INFO"))

    raw_dir = Path(cfg["paths"]["data_raw"])
    raw_dir.mkdir(parents=True, exist_ok=True)

    gdc_cfg = cfg.get("gdc", {})
    project = gdc_cfg.get("project_id") or cfg.get("project", "TCGA-LUAD")
    conda_prefix = Path(gdc_cfg.get("conda_prefix", "~/mnt/datapool/conda-envs/tcgabiolinks")).expanduser().resolve()

    rscript = conda_prefix / "bin" / "Rscript"
    if not rscript.exists():
        raise FileNotFoundError(
            f"TCGAbiolinks environment not found at {rscript}. "
            f"Set gdc.conda_prefix in config.yaml or create the environment first."
        )

    env = os.environ.copy()
    env["PATH"] = f"{conda_prefix / 'bin'}:{env.get('PATH', '')}"
    env["GDC_RSCRIPT"] = str(rscript)
    # Persist the known-working sesame workaround unless the caller overrides it.
    env.setdefault("R_LIBS_USER", "/tmp/tcgabiolinks-rlib")

    clinical_path = raw_dir / f"{project}.clinical.tsv.gz"
    survival_path = raw_dir / f"{project}.survival.tsv.gz"
    rnaseq_path = raw_dir / f"{project}.htseq_counts.tsv.gz"
    mutation_path = raw_dir / f"{project}.mutect2_snv.tsv.gz"
    cnv_path = raw_dir / f"{project}.gistic.tsv.gz"
    methylation_path = raw_dir / f"{project}.methylation450.tsv.gz"
    mirna_path = raw_dir / f"{project}.mirna.tsv.gz"

    script = textwrap.dedent(
        f'''
        suppressPackageStartupMessages({{
          library(TCGAbiolinks)
          library(SummarizedExperiment)
        }})

        project <- "{project}"
        raw_dir <- "{str(raw_dir)}"
        raw_real <- normalizePath(raw_dir, winslash = "/", mustWork = FALSE)
        dir.create(raw_real, recursive = TRUE, showWarnings = FALSE)
        setwd(raw_real)
        options(tmpdir = raw_real)

        write_gz_tsv <- function(df, path) {{
          df <- as.data.frame(df, stringsAsFactors = FALSE)
          for (nm in names(df)) {{
            if (is.list(df[[nm]])) {{
              df[[nm]] <- vapply(df[[nm]], function(x) paste(as.character(x), collapse = ";"), character(1))
            }}
          }}
          con <- gzfile(path, "wt")
          on.exit(close(con), add = TRUE)
          write.table(df, con, sep = "\t", quote = FALSE, row.names = FALSE)
        }}

        save_object <- function(obj, path, row_id_name="feature", label="object") {{
          message(sprintf("CLASS[%s]: %s", label, paste(class(obj), collapse=",")))
          if (inherits(obj, "SummarizedExperiment")) {{
            mat <- assay(obj)
            df <- data.frame(row_id = rownames(mat), mat, check.names = FALSE)
            names(df)[1] <- row_id_name
          }} else if (is.matrix(obj)) {{
            df <- data.frame(row_id = rownames(obj), obj, check.names = FALSE)
            names(df)[1] <- row_id_name
          }} else if (is.data.frame(obj)) {{
            df <- obj
          }} else {{
            stop(sprintf("Unsupported object class for %s: %s", label, paste(class(obj), collapse=",")))
          }}
          for (nm in names(df)) {{
            if (is.list(df[[nm]])) {{
              df[[nm]] <- vapply(df[[nm]], function(x) paste(as.character(x), collapse = ";"), character(1))
            }}
          }}
          con <- gzfile(path, "wt")
          on.exit(close(con), add = TRUE)
          write.table(df, con, sep = "\t", quote = FALSE, row.names = FALSE)
        }}

        file_ready <- function(path) {{
          file.exists(path) && file.info(path)$size > 0
        }}

        # Clinical + survival
        get_clinical_safe <- function(project) {{
          tryCatch({{
            GDCquery_clinic(project = project, type = "clinical")
          }}, error = function(e) {{
            message(sprintf("WARN[clinical]: GDCquery_clinic failed (%s); falling back to indexed clinical query", e$message))
            q_clin <- GDCquery(
              project = project,
              data.category = "Clinical",
              data.type = "Clinical Supplement",
              data.format = "BCR Biotab"
            )
            GDCdownload(q_clin, method = "api", files.per.chunk = 20, directory = raw_real)
            clin_files <- list.files(raw_real, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
            clin_files <- clin_files[grepl("patient", basename(clin_files), ignore.case = TRUE) & grepl("clinical|nationwidechildrens", basename(clin_files), ignore.case = TRUE)]
            if (length(clin_files) == 0) {{
              stop(sprintf("No downloaded clinical patient file found for fallback clinical query under %s", raw_real))
            }}
            clin <- read.delim(clin_files[1], sep = "\t", quote = "", comment.char = "", check.names = FALSE, stringsAsFactors = FALSE)
            as.data.frame(clin, stringsAsFactors = FALSE)
          }})
        }}

        if (!file_ready("{clinical_path}")) {{
          clin <- get_clinical_safe(project)
          write_gz_tsv(clin, "{clinical_path}")
        }} else {{
          message("SKIP[clinical]: output exists")
        }}

        if (!file_ready("{survival_path}")) {{
          surv <- tryCatch({{
            GDCquery_clinic(project = project, type = "survival")
          }}, error = function(e) NULL)
          if (!is.null(surv)) {{
            write_gz_tsv(surv, "{survival_path}")
          }}
        }} else {{
          message("SKIP[survival]: output exists")
        }}

        safe_download_prepare <- function(query, label, summarized = TRUE) {{
          message(sprintf("WORKDIR[%s]: %s", label, getwd()))
          if (is.null(query$results[[1]]) || nrow(query$results[[1]]) == 0) {{
            message(sprintf("SKIP[%s]: query returned no results", label))
            return(NULL)
          }}
          GDCdownload(query, method = "api", files.per.chunk = 20, directory = raw_real)
          tryCatch({{
            GDCprepare(query, directory = raw_real, summarizedExperiment = summarized)
          }}, error = function(e) {{
            msg <- conditionMessage(e)
            if (grepl("disease_response", msg, fixed = TRUE)) {{
              message(sprintf("WARN[%s]: retrying GDCprepare with stripped clinical metadata", label))
              old_colDataPrepare <- get("colDataPrepare", envir = asNamespace("TCGAbiolinks"))
              assign("colDataPrepare", function(cases) {{
                cd <- old_colDataPrepare(cases)
                if ("disease_response" %in% colnames(cd)) {{
                  cd$disease_response <- NULL
                }}
                cd
              }}, envir = asNamespace("TCGAbiolinks"))
              on.exit(assign("colDataPrepare", old_colDataPrepare, envir = asNamespace("TCGAbiolinks")), add = TRUE)
              return(GDCprepare(query, directory = raw_real, summarizedExperiment = summarized))
            }}
            stop(e)
          }})
        }}

        dedup_query_cases <- function(query, label) {{
          res <- query$results[[1]]
          if (is.null(res) || nrow(res) == 0) {{
            return(query)
          }}
          if ("sample.submitter_id" %in% colnames(res)) {{
            multi_sample <- grepl(";", res$sample.submitter_id, fixed = TRUE)
            if (any(multi_sample, na.rm = TRUE)) {{
              message(sprintf("DEDUP[%s]: dropping %d multi-sample rows", label, sum(multi_sample, na.rm = TRUE)))
              res <- res[!multi_sample, , drop = FALSE]
            }}
          }}
          if ("sample_type" %in% colnames(res)) {{
            primary_idx <- res$sample_type == "Primary Tumor"
            primary_idx[is.na(primary_idx)] <- FALSE
            if (any(primary_idx)) {{
              message(sprintf("DEDUP[%s]: keeping %d Primary Tumor rows", label, sum(primary_idx)))
              res <- res[primary_idx, , drop = FALSE]
            }}
          }}
          if ("sample.submitter_id" %in% colnames(res)) {{
            dup_sample <- duplicated(res$sample.submitter_id)
            if (any(dup_sample, na.rm = TRUE)) {{
              message(sprintf("DEDUP[%s]: dropping %d duplicate sample rows", label, sum(dup_sample, na.rm = TRUE)))
              res <- res[!dup_sample, , drop = FALSE]
            }}
          }}
          if (nrow(res) == 0) {{
            message(sprintf("SKIP[%s]: no rows remain after dedup", label))
            query$results[[1]] <- res
            return(query)
          }}
          query$results[[1]] <- res
          query
        }}

        # RNA-seq counts
        if (!file_ready("{rnaseq_path}")) {{
          q_rna <- GDCquery(
            project = project,
            data.category = "Transcriptome Profiling",
            data.type = "Gene Expression Quantification",
            workflow.type = "STAR - Counts"
          )
          se_rna <- safe_download_prepare(q_rna, "rna")
          save_object(se_rna, "{rnaseq_path}", row_id_name = "Ensembl_ID", label = "rna")
        }} else {{
          message("SKIP[rna]: output exists")
        }}

        # miRNA
        if (!file_ready("{mirna_path}")) {{
          q_mir <- GDCquery(
            project = project,
            data.category = "Transcriptome Profiling",
            data.type = "miRNA Expression Quantification"
          )
          se_mir <- safe_download_prepare(q_mir, "mirna")
          save_object(se_mir, "{mirna_path}", row_id_name = "miRNA_ID", label = "mirna")
        }} else {{
          message("SKIP[mirna]: output exists")
        }}

        # Methylation 450K
        if (!file_ready("{methylation_path}")) {{
          q_meth <- GDCquery(
            project = project,
            data.category = "DNA Methylation",
            data.type = "Methylation Beta Value",
            platform = "Illumina Human Methylation 450"
          )
          se_meth <- safe_download_prepare(q_meth, "methylation")
          save_object(se_meth, "{methylation_path}", row_id_name = "Probe", label = "methylation")
        }} else {{
          message("SKIP[methylation]: output exists")
        }}

        # CNV gene-level copy number
        if (!file_ready("{cnv_path}")) {{
          q_cnv <- GDCquery(
            project = project,
            data.category = "Copy Number Variation",
            data.type = "Gene Level Copy Number"
          )
          q_cnv <- dedup_query_cases(q_cnv, "cnv")
          cnv_obj <- safe_download_prepare(q_cnv, "cnv", summarized = FALSE)
          if (!is.null(cnv_obj)) {{
            save_object(cnv_obj, "{cnv_path}", row_id_name = "Gene", label = "cnv")
          }}
        }} else {{
          message("SKIP[cnv]: output exists")
        }}

        # Mutations (current GDC SNV workflow)
        if (!file_ready("{mutation_path}")) {{
          q_mut <- GDCquery(
            project = project,
            data.category = "Simple Nucleotide Variation",
            data.type = "Masked Somatic Mutation",
            workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
          )
          mut_obj <- safe_download_prepare(q_mut, "mutation", summarized = FALSE)
          if (!is.null(mut_obj)) {{
            save_object(mut_obj, "{mutation_path}", row_id_name = "gene", label = "mutation")
          }}
        }} else {{
          message("SKIP[mutation]: output exists")
        }}
        '''
    )

    log.info(f"Using TCGAbiolinks from {conda_prefix}")
    _run_r(script, env=env, log=log)

    out = {}
    for name, path in {
        "clinical": clinical_path,
        "survival": survival_path,
        "rnaseq": rnaseq_path,
        "mutation": mutation_path,
        "cnv": cnv_path,
        "methylation": methylation_path,
        "mirna": mirna_path,
    }.items():
        if path.exists() and path.stat().st_size > 0:
            out[name] = path

    log.info(f"Downloaded / prepared {len(out)} modalities into {raw_dir}")
    return out
