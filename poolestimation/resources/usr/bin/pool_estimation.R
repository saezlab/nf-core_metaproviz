#!/usr/bin/env Rscript
# =============================================================================
# pool_estimation.R
#
# Run MetaProViz::pool_estimation on a SummarizedExperiment to assess metabolite
# detection quality from pool samples (mixtures of all experimental samples
# measured several times during the LC-MS run).
#
# Pool samples can be identified in three ways (top precedence first):
#
#   1) Explicit list           --pool_samples "POOL1,POOL2,..."
#   2) Metadata-column lookup  --pool_metadata_col Conditions
#                              --pool_metadata_value Pool
#   3) Pattern match on row    --pool_pattern "pool"   (default; case-insensitive
#      names                     substring match on sample IDs)
#
# Pool-estimation parameters:
#   --cutoff_cv          CV % above which a feature is flagged HighVar (default: 30).
#   --filter_high_var    TRUE  >>> drop HighVar features from the SE (default).
#                        FALSE >>> keep them; only report.
#
# Outputs (all named <prefix>.<suffix>; --prefix defaults to "pool_estimation"):
#   <prefix>.pool_filtered.rds  SE with HighVar features dropped (or copy)
#   <prefix>.cv.tsv             full CV table
#   <prefix>.high_var.txt       one feature name per line
#   <prefix>.plots.rds          list of ggplot objects (so you can inspect /
#                               reuse them in R: readRDS("...")$PCAPlot…)
#   <prefix>.report.html        self-contained HTML report (PNGs embedded
#                               base64; no extra files needed)
#   <prefix>.log                INFO/WARN log
# =============================================================================

set.seed(42)
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(magrittr)
  library(ggplot2)
  library(SummarizedExperiment)
  library(S4Vectors)
})

suppressPackageStartupMessages(library(MetaProViz))

# ─────────────────────────────────────────────────────────────────────────────
# ARGUMENT PARSING
# ─────────────────────────────────────────────────────────────────────────────

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  val <- args[idx + 1]
  if (startsWith(val, "--")) return(default)
  val
}

se_path             <- get_arg("--se")
data_matrix_path    <- get_arg("--data_matrix")
feature_matrix_path <- get_arg("--feature_matrix")
sample_matrix_path  <- get_arg("--sample_matrix")
pool_pattern        <- get_arg("--pool_pattern", "pool")
pool_samples_str    <- get_arg("--pool_samples", "")
pool_metadata_col   <- get_arg("--pool_metadata_col", "")
pool_metadata_value <- get_arg("--pool_metadata_value", "")
cutoff_cv           <- as.numeric(get_arg("--cutoff_cv", "30"))
filter_high_var_str <- toupper(get_arg("--filter_high_var", "TRUE"))
filter_high_var     <- filter_high_var_str == "TRUE"
# All output filenames are `<prefix>.<suffix>`, per nf-core naming
# convention (https://nf-co.re/docs/specifications/components/modules/naming-conventions):
# "output file names SHOULD consist of only ${prefix} and the file-format
# suffix" — needed so filenames don't collide across samples when this
# module runs on many samples in one pipeline.
prefix              <- get_arg("--prefix", "pool_estimation")
se_out              <- get_arg("--se_out", paste0(prefix, ".pool_filtered.rds"))
cv_out              <- paste0(prefix, ".cv.tsv")
high_var_out        <- paste0(prefix, ".high_var.txt")
plots_out           <- paste0(prefix, ".plots.rds")
report_out          <- paste0(prefix, ".report.html")
log_out             <- paste0(prefix, ".log")

# ─────────────────────────────────────────────────────────────────────────────
# LOGGING
# ─────────────────────────────────────────────────────────────────────────────

log_lines  <- character(0)
warn_count <- 0; err_count <- 0

log_msg <- function(level = "INFO", ...) {
  txt <- paste0("[", level, "] ", paste(..., sep = ""))
  log_lines <<- c(log_lines, txt); message(txt)
  if (level == "WARN")  warn_count <<- warn_count + 1
  if (level == "ERROR") err_count  <<- err_count  + 1
}
log_section <- function(title) {
  sep <- paste(rep("─", 60), collapse = "")
  log_lines <<- c(log_lines, "", sep, paste0("  ", title), sep)
  message(sep, "\n  ", title, "\n", sep)
}
flush_log <- function(path = log_out) {
  writeLines(c(log_lines, "",
               sprintf("SUMMARY: %d warning(s), %d error(s)",
                       warn_count, err_count)), path)
}
abort <- function(...) { log_msg("ERROR", ...); flush_log(); stop(paste(...), call. = FALSE) }

# ─────────────────────────────────────────────────────────────────────────────
# LOAD SE
# ─────────────────────────────────────────────────────────────────────────────

log_section("Loading input data")

# Two acceptable input shapes, mutually exclusive: --se (a SummarizedExperiment
# .rds), or all three of --data_matrix/--feature_matrix/--sample_matrix (plain
# TSV/CSV, language-independent). Test data is provided in both shapes; this
# module accepts either so users aren't forced into one or the other.
have_se  <- !is.null(se_path) && nzchar(se_path)
have_csv <- !is.null(data_matrix_path) && nzchar(data_matrix_path) &&
            !is.null(feature_matrix_path) && nzchar(feature_matrix_path) &&
            !is.null(sample_matrix_path) && nzchar(sample_matrix_path)
have_partial_csv <- !have_csv && (
  (!is.null(data_matrix_path) && nzchar(data_matrix_path)) ||
  (!is.null(feature_matrix_path) && nzchar(feature_matrix_path)) ||
  (!is.null(sample_matrix_path) && nzchar(sample_matrix_path))
)

if (have_se && have_csv)
  abort("Provide either --se, or all three of --data_matrix/--feature_matrix/",
        "--sample_matrix, not both.")
if (have_partial_csv)
  abort("--data_matrix, --feature_matrix, and --sample_matrix must all be ",
        "provided together.")
if (!have_se && !have_csv)
  abort("Provide either --se <path>, or all three of --data_matrix/",
        "--feature_matrix/--sample_matrix.")

if (have_se) {
  if (!file.exists(se_path)) abort("SE file not found: ", se_path)
  se <- tryCatch(readRDS(se_path),
                 error = function(e) abort("Failed to read SE: ",
                                           conditionMessage(e)))
  if (!is(se, "SummarizedExperiment"))
    abort("Object in ", se_path, " is not a SummarizedExperiment.")
  log_msg("INFO", "Loaded SE from ", basename(se_path), ": ", nrow(se),
          " features × ", ncol(se), " samples.")
} else {
  for (p in c(data_matrix_path, feature_matrix_path, sample_matrix_path))
    if (!file.exists(p)) abort("Input file not found: ", p)

  read_flat <- function(path) {
    sep <- if (grepl("\\.csv$", path, ignore.case = TRUE)) "," else "\t"
    tryCatch(read.delim(path, sep = sep, check.names = FALSE,
                        stringsAsFactors = FALSE),
             error = function(e) abort("Failed to read ", path, ": ",
                                       conditionMessage(e)))
  }

  dm <- read_flat(data_matrix_path)
  fm <- read_flat(feature_matrix_path)
  sm <- read_flat(sample_matrix_path)

  if (ncol(dm) < 2)
    abort("--data_matrix must have a feature-ID column plus at least one ",
          "sample column.")

  feature_ids <- as.character(dm[[1]])
  assay_matrix <- as.matrix(dm[, -1, drop = FALSE])
  rownames(assay_matrix) <- feature_ids
  mode(assay_matrix) <- "numeric"

  sample_ids <- as.character(sm[[1]])
  col_data <- sm[, -1, drop = FALSE]
  rownames(col_data) <- sample_ids

  feature_ids2 <- as.character(fm[[1]])
  row_data <- fm[, -1, drop = FALSE]
  rownames(row_data) <- feature_ids2

  missing_samples <- setdiff(colnames(assay_matrix), sample_ids)
  if (length(missing_samples) > 0)
    abort("--data_matrix has sample column(s) not present in ",
          "--sample_matrix: ", paste(missing_samples, collapse = ", "))
  missing_features <- setdiff(rownames(assay_matrix), feature_ids2)
  if (length(missing_features) > 0)
    abort("--data_matrix has feature(s) not present in --feature_matrix: ",
          paste(head(missing_features, 5), collapse = ", "),
          if (length(missing_features) > 5) ", ..." else "")

  col_data <- col_data[colnames(assay_matrix), , drop = FALSE]
  row_data <- row_data[rownames(assay_matrix), , drop = FALSE]

  se <- SummarizedExperiment(assays = list(counts = assay_matrix),
                             colData = col_data, rowData = row_data)
  log_msg("INFO", "Built SE from data_matrix/feature_matrix/sample_matrix: ",
          nrow(se), " features × ", ncol(se), " samples.")
}

# Used in the report/log wherever the input source needs naming; se_path is
# NULL when the module was given the 3 flat files instead of an .rds.
source_label <- if (have_se) {
  basename(se_path)
} else {
  paste(basename(data_matrix_path), basename(feature_matrix_path),
        basename(sample_matrix_path), sep = " + ")
}

# Data layout for MetaProViz: samples x features (rownames = sample IDs)
assay_df <- t(assay(se, 1)) %>% as.data.frame(check.names = FALSE)
sample_info_df <- as.data.frame(colData(se), check.names = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# IDENTIFY POOL SAMPLES
# ─────────────────────────────────────────────────────────────────────────────

log_section("Identifying pool samples")

pool_ids <- character(0)
strategy <- "none"

# 1) Explicit list ────────────────────────────────────────────────────────────
if (nzchar(pool_samples_str)) {
  asked <- trimws(strsplit(pool_samples_str, ",")[[1]])
  asked <- asked[nzchar(asked)]
  pool_ids <- intersect(asked, rownames(assay_df))
  missing  <- setdiff(asked, rownames(assay_df))
  if (length(missing) > 0)
    log_msg("WARN", "Pool sample(s) not found in SE: ",
            paste(missing, collapse = ", "))
  strategy <- "explicit_list"
}

# 2) Metadata-column lookup ──────────────────────────────────────────────────
if (length(pool_ids) == 0 && nzchar(pool_metadata_col)) {
  if (!pool_metadata_col %in% colnames(sample_info_df)) {
    log_msg("WARN", "pool_metadata_col '", pool_metadata_col,
            "' not in colData. Available: ",
            paste(colnames(sample_info_df), collapse = ", "))
  } else if (!nzchar(pool_metadata_value)) {
    log_msg("WARN", "pool_metadata_col was set without --pool_metadata_value. ",
            "Skipping metadata-column strategy.")
  } else {
    hit <- which(as.character(sample_info_df[[pool_metadata_col]]) ==
                 pool_metadata_value)
    pool_ids <- rownames(sample_info_df)[hit]
    strategy <- "metadata_column"
  }
}

# 3) Pattern on rownames (default fallback) ──────────────────────────────────
if (length(pool_ids) == 0 && nzchar(pool_pattern)) {
  hit <- grepl(pool_pattern, rownames(assay_df), ignore.case = TRUE)
  pool_ids <- rownames(assay_df)[hit]
  strategy <- paste0("pattern '", pool_pattern, "'")
}

if (length(pool_ids) == 0) {
  abort("No pool samples found. Either set --pool_samples explicitly, ",
        "or --pool_metadata_col + --pool_metadata_value, or ensure your ",
        "pool sample IDs match the pattern (default: 'pool').")
}

log_msg("INFO", "Pool detection strategy: ", strategy)
log_msg("INFO", "Pool samples (", length(pool_ids), "): ",
        paste(head(pool_ids, 12), collapse = ", "),
        if (length(pool_ids) > 12) " ..." else "")

# ─────────────────────────────────────────────────────────────────────────────
# RUN POOL ESTIMATION
# ─────────────────────────────────────────────────────────────────────────────

log_section("Running MetaProViz::pool_estimation")

# Some MetaProViz versions require a writable log dir; create a per-run scratch.
tmp_log_dir <- tempfile("metaproviz_log_")
dir.create(tmp_log_dir, recursive = TRUE, showWarnings = FALSE)
old_wd <- getwd(); on.exit(setwd(old_wd), add = TRUE)

# MetaProViz expects:
#   - a Conditions column in metadata_sample
#   - PoolSamples = the value inside that column that marks pool rows
#     (not a separate column name).
#
# To make this robust regardless of what's in the input:
#   - If colData has no Conditions column, create one (all "Sample" by default).
#   - Save the original Conditions to Conditions_original so nothing is lost.
#   - Stamp the pool sample rows with the literal "Pool" marker.
#   - Warn loudly if a non-pool sample already had "Pool" in Conditions, since
#     that would silently make MetaProViz treat it as a pool.

pool_label <- "Pool"

if (!"Conditions" %in% colnames(sample_info_df)) {
  log_msg("INFO", "No 'Conditions' column found in colData — creating one.")
  sample_info_df$Conditions <- "Sample"
} else {
  sample_info_df$Conditions_original <- sample_info_df$Conditions
}

is_pool_row    <- rownames(sample_info_df) %in% pool_ids
already_pool   <- sample_info_df$Conditions == pool_label & !is_pool_row
if (any(already_pool, na.rm = TRUE)) {
  log_msg("WARN",
    sum(already_pool), " non-pool sample(s) already had '", pool_label,
    "' in their Conditions column: ",
    paste(head(rownames(sample_info_df)[already_pool], 5), collapse = ", "),
    if (sum(already_pool) > 5) " ..." else "",
    ". MetaProViz would treat them as pools — overwriting their Conditions ",
    "with 'Sample' to avoid this. Original kept in Conditions_original.")
  sample_info_df$Conditions[already_pool] <- "Sample"
}

sample_info_df$Conditions[is_pool_row] <- pool_label
log_msg("INFO", "Stamped Conditions = '", pool_label, "' on ",
        sum(is_pool_row), " pool sample row(s).")

# Sanity log: show distribution of Conditions before calling MetaProViz.
log_msg("INFO", "Conditions levels going into pool_estimation: ",
        paste(sprintf("%s=%d", names(table(sample_info_df$Conditions)),
                      as.integer(table(sample_info_df$Conditions))),
              collapse = ", "))

pe <- tryCatch(
  MetaProViz::pool_estimation(
    data            = assay_df,
    metadata_sample = sample_info_df,
    metadata_info   = c(PoolSamples = pool_label, Conditions = "Conditions"),
    cutoff_cv       = cutoff_cv,
    print_plot      = FALSE,
    save_plot       = NULL
  ),
  error = function(e) abort("pool_estimation() failed: ", conditionMessage(e))
)

# ─────────────────────────────────────────────────────────────────────────────
# EXTRACT RESULTS
# ─────────────────────────────────────────────────────────────────────────────

log_section("Extracting results")

cv_df <- pe[["DF"]][["CV"]]
if (is.null(cv_df))
  abort("pool_estimation() did not return DF$CV — unexpected result shape.")
log_msg("INFO", "CV table: ", nrow(cv_df), " feature(s), columns: ",
        paste(colnames(cv_df), collapse = ", "))

high_var <- cv_df %>% filter(HighVar == TRUE) %>% pull(Metabolite)
log_msg("INFO", "High-variance features at CV > ", cutoff_cv,
        "%: ", length(high_var), " / ", nrow(cv_df))

plots <- pe[["Plot"]]
if (is.null(plots) || length(plots) == 0) {
  log_msg("WARN", "pool_estimation() did not return any plots.")
} else {
  log_msg("INFO", "Plots returned: ", paste(names(plots), collapse = ", "))
}

# ─────────────────────────────────────────────────────────────────────────────
# OPTIONAL FILTERING
# ─────────────────────────────────────────────────────────────────────────────

log_section("Filtering high-variance features (optional)")

if (filter_high_var && length(high_var) > 0) {
  keep <- setdiff(rownames(se), high_var)
  log_msg("INFO", "Removing ", length(high_var),
          " high-variance feature(s); keeping ", length(keep), ".")
  se_filtered <- se[keep, ]
  metadata(se_filtered)$pool_estimation <- list(
    cutoff_cv      = cutoff_cv,
    pool_samples   = pool_ids,
    high_var_feats = high_var,
    strategy       = strategy
  )
} else {
  se_filtered <- se
  if (filter_high_var) {
    log_msg("INFO", "No HighVar features to remove — SE unchanged.")
  } else {
    log_msg("INFO", "filter_high_var=FALSE — SE unchanged.")
  }
  metadata(se_filtered)$pool_estimation <- list(
    cutoff_cv    = cutoff_cv,
    pool_samples = pool_ids,
    strategy     = strategy
  )
}

# ─────────────────────────────────────────────────────────────────────────────
# WRITE OUTPUTS
# ─────────────────────────────────────────────────────────────────────────────

log_section("Writing outputs")

saveRDS(se_filtered, se_out);                       log_msg("INFO", "Written: ", se_out)

write.table(cv_df, cv_out,
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
log_msg("INFO", "Written: ", cv_out)

writeLines(if (length(high_var) > 0) high_var else "",
           high_var_out)
log_msg("INFO", "Written: ", high_var_out)

saveRDS(plots, plots_out)
log_msg("INFO", "Written: ", plots_out, " (a named list of ggplots)")

# ─────────────────────────────────────────────────────────────────────────────
# HTML REPORT
# Self-contained: ggplots >>> PNG >>> base64 >>> embedded <img>.
# ─────────────────────────────────────────────────────────────────────────────

log_section("Building HTML report")

html_escape <- function(x) {
  x <- gsub("&", "&amp;",  x, fixed = TRUE)
  x <- gsub("<", "&lt;",   x, fixed = TRUE)
  x <- gsub(">", "&gt;",   x, fixed = TRUE)
  x
}
plot_to_b64 <- function(p, width = 8, height = 5, dpi = 100) {
  tmp <- tempfile(fileext = ".png")
  ggplot2::ggsave(tmp, plot = p, width = width, height = height,
                  dpi = dpi, units = "in", bg = "white")
  raw <- readBin(tmp, what = "raw", n = file.info(tmp)$size)
  unlink(tmp)
  paste0("data:image/png;base64,",
         base64enc::base64encode(raw))
}
df_to_html_table <- function(df, max_rows = 5, source_file = NULL) {
  # Big tables now show only the top `max_rows` rows; the full data lives in
  # the TSV next to the HTML. Tables are also wrapped in a scrollable
  # container with a sticky header for any inline scrolling.
  rows_shown <- head(df, max_rows)
  th <- paste0("<th>", html_escape(colnames(rows_shown)), "</th>", collapse = "")
  body_rows <- apply(rows_shown, 1, function(r) {
    paste0("<tr>",
           paste0("<td>", html_escape(as.character(r)), "</td>", collapse = ""),
           "</tr>")
  })
  more <- if (nrow(df) > max_rows) {
    src_hint <- if (!is.null(source_file))
      sprintf(" (full data in %s)", html_escape(source_file)) else ""
    sprintf("<p class='more-rows'><em>Showing first %d of %d row(s)%s.</em></p>",
            max_rows, nrow(df), src_hint)
  } else {
    ""
  }
  paste0(
    '<div class="scroll-table"><table><thead><tr>', th,
    "</tr></thead><tbody>",
    paste(body_rows, collapse = ""), "</tbody></table></div>", more
  )
}

# Build embedded plot section
plot_html <- ""
if (!is.null(plots) && length(plots) > 0) {
  for (nm in names(plots)) {
    p <- plots[[nm]]
    if (!inherits(p, "ggplot")) next
    src <- tryCatch(plot_to_b64(p),
                    error = function(e) {
                      log_msg("WARN", "Could not render plot '", nm, "': ",
                              conditionMessage(e)); NA_character_
                    })
    if (!is.na(src)) {
      plot_html <- paste0(
        plot_html,
        sprintf('<h3>%s</h3>\n<img src="%s" alt="%s" />\n',
                html_escape(nm), src, html_escape(nm))
      )
    }
  }
}

high_var_html <- if (length(high_var) > 0) {
  paste0("<ul>",
         paste0("<li>", html_escape(high_var), "</li>", collapse = ""),
         "</ul>")
} else "<p><em>None — all features below the CV cutoff.</em></p>"

run_meta_html <- sprintf(
  "<dl>
     <dt>Input source</dt><dd>%s</dd>
     <dt>n features (input)</dt><dd>%d</dd>
     <dt>n samples (input)</dt><dd>%d</dd>
     <dt>Pool detection strategy</dt><dd>%s</dd>
     <dt>Pool samples (%d)</dt><dd>%s</dd>
     <dt>CV cutoff</dt><dd>%s%%</dd>
     <dt>HighVar features</dt><dd>%d / %d</dd>
     <dt>Filter HighVar from SE</dt><dd>%s</dd>
   </dl>",
  html_escape(source_label),
  nrow(se), ncol(se),
  html_escape(strategy),
  length(pool_ids),
  html_escape(paste0(
    paste(head(pool_ids, 25), collapse = ", "),
    if (length(pool_ids) > 25) sprintf(" (+%d more)", length(pool_ids) - 25) else ""
  )),
  format(cutoff_cv),
  length(high_var), nrow(cv_df),
  if (filter_high_var) "TRUE" else "FALSE"
)

html <- sprintf('<!DOCTYPE html>
<html><head><meta charset="utf-8"><title>Pool Estimation Report — %s</title>
<style>
  body  { font-family: Arial, sans-serif; max-width: 980px; margin: 2em auto;
          color: #1e1e1e; line-height: 1.5; padding: 0 1em; }
  h1    { color: #1F3A5F; border-bottom: 2px solid #1F3A5F; padding-bottom: .3em; }
  h2    { color: #1F3A5F; margin-top: 2em; }
  h3    { color: #345; margin-top: 1.4em; }
  dl    { display: grid; grid-template-columns: max-content 1fr; gap: .3em 1em; }
  dt    { font-weight: bold; color: #345; }
  table { border-collapse: collapse; width: 100%%; font-size: .9em; }
  th,td { border: 1px solid #ccc; padding: .35em .6em; text-align: left; }
  th    { background: #f0f4f8; }
  tr:nth-child(even) td { background: #fafafa; }
  img   { max-width: 100%%; border: 1px solid #ddd; padding: 4px;
          background: #fff; margin: .4em 0 1.2em; }
  code  { background: #f2f2f2; padding: .15em .35em; border-radius: 3px; }
    pre   { background: #f7f9fc; border: 1px solid #d9e2ec; border-radius: 6px;
      padding: .7em .9em; overflow-x: auto; }
    li    { margin: .25em 0; }
  ul    { margin-top: .3em; }
  /* Scrollable table container with sticky header */
  .scroll-table       { max-height: 320px; overflow-y: auto;
                        border: 1px solid #cbd6e4; border-radius: 6px;
                        margin: .4em 0 .6em; }
  .scroll-table table { border: 0; margin: 0; }
  .scroll-table thead th { position: sticky; top: 0; z-index: 1;
                           background: #f0f4f8; }
  .more-rows { color: #555; font-size: .9em; margin: .2em 0 1.2em; }
</style></head><body>

<h1>Pool Estimation Report</h1>
<p><em>Generated by the nf_metabolism pipeline (POOL_ESTIMATION module).</em></p>

<h2>Run summary</h2>
%s

<h2>Citations</h2>
<p>
If you use results or plots from this analysis in a publication, please cite:
</p>
<ul>
  <li><strong>MetaProViz:</strong> Please cite:
  <a href="https://doi.org/10.1038/s44320-026-00231-8">Schmidt et al., Integrated metabolomics data analysis to generate mechanistic hypotheses with MetaProViz, Molecular Systems Biology 2026.</a>.</li>
  <li><strong>Dependencies:</strong> None for this module.</li>
</ul>

<h2>Description</h2>
<p>
This report uses pooled quality-control (QC) samples to assess technical data
quality. Pool samples are a homogeneous mixture of all study samples, so they
should be stable across repeated injections. The result of
<code>Pool_Estimation()</code> is a feature-level CV table; features with high
variability can be considered for removal from downstream analyses.
</p>
<p>
QC-based identification of low-quality features is multifaceted, and pool-sample
QC is one common approach. Specifically, the coefficient of variation (CV) helps
identify features with high variability and therefore low detection accuracy.
As a practical rule of thumb, pool-sample CV values are typically expected to be
below <strong>%s%%</strong> (default cutoff).
</p>
<p>
CV is calculated as:
</p>
<p><code>CV(%%) = 100 × sd(x1, x2, ..., xn) / mean(x1, x2, ..., xn)</code></p>
<p>
This module computes CVs through MetaProViz and extracts the returned CV table.
</p>

<h2>Plots</h2>
<p>
<strong>How to read the plots:</strong>
</p>
<ul>
  <li><strong>PCA plot:</strong> Pool samples are expected to cluster near the coordinate origin,
  reflecting their role as a homogeneous mixture of all samples. This pattern is
  most interpretable when all sample groups are included in the PCA view.</li>
  <li><strong>Histogram:</strong> x-axis = CV and y-axis = frequency. The cutoff line is set by
  the user (typically %s%%). Most features are expected at lower CV values (left side),
  with frequencies tapering toward higher CV values.</li>
  <li><strong>Violin plot:</strong> Distribution of metabolite CV values split into two
  groups relative to the cutoff (%s%%): features below cutoff vs features above
  cutoff (high-variance candidates).</li>
</ul>
<p>
Use these diagnostics together with the high-variance feature list and full CV
table to decide whether high-CV features should be removed.
</p>
%s

<h2>High-variance features (CV > %s%%)</h2>
%s

<h2>Full CV table</h2>
%s

<h2>Run notes (log)</h2>
<p>A separate file, <code>%s</code>, records every step of this run,
including any warnings. Check it first if a result here looks unexpected.</p>

</body></html>',
  html_escape(source_label),
  run_meta_html,
  format(cutoff_cv),
  format(cutoff_cv),
  format(cutoff_cv),
  if (nzchar(plot_html)) plot_html else "<p><em>No plots returned.</em></p>",
  format(cutoff_cv),
  high_var_html,
  # Show the FULL CV table in HTML — it's wrapped in a scrollable container
  # (.scroll-table) so it stays compact. The PDF build collapses it to the
  # header + a "see full data" note (handled by combine_reports.py).
  df_to_html_table(cv_df, max_rows = 5000,
                   source_file = cv_out),
  html_escape(log_out)
)
writeLines(html, report_out)
log_msg("INFO", "Written: ", report_out)

flush_log(log_out)
log_msg("INFO", "Written: ", log_out)

message(sprintf(
  "\n✓ pool_estimation complete — %d / %d features flagged HighVar | %d warning(s) | %d error(s)",
  length(high_var), nrow(cv_df), warn_count, err_count
))
