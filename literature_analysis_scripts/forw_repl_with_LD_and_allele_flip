# ==============================================================================
# FORWARD LITERATURE REVIEW PIPELINE (copy-paste ready)
# - 5 replication classes
# - LD mode: abs / signed_raw / signed_corrected
# - LD threshold rule: abs (|LD|>=thr) / signed (LD>=thr)
# - Distinguishes LD missing due to variants absent in preQC (.bim)
#   vs LD absent from LD file (filtered because too low) while variants exist in preQC
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(scales)
  library(data.table)
})

# ==============================================================================
# SETTINGS
# ==============================================================================
thr       <- -log10(5e-8)
LD_THR    <- 0.8             # minimum LD for replication
REGION_BP <- 1e6
ANCHOR_COL <- "SNPID_LEAD"

# Forward literature review input location
base_dir   <- "/exchange/healthds/pQTL/BELIEVE/gwasstudio_output"

# BELIEVE preQC genotype bfiles directory (chr*.bim used to check variant presence)
bfile_dir  <- "/scratch/laura.bondi/BELIEVE_bfile"

# BELIEVE LD files (PLINK --r output)
ld_base_dir <- "/scratch/laura.bondi/PLINK_LD_RESUME_FINAL"
ld_file_for_chr <- function(chr) file.path(ld_base_dir, paste0("chr", chr), paste0("chr", chr, "_r.ld.gz"))

# Column names in the per-study CSVs
LEAD_ID_COL  <- "SNPID_LEAD"
EXACT_ID_COL <- "SNPID_EXACT"

# --- CHOICES ON HOW TO DEFINE REPLICATION ---
# I choose use LD without sign 
# If using LD with sign, one must be careful about sign changes do to allele flip and use signed_corrected
# LD sign handling:
#   "abs"             -> use abs(r)
#   "signed_raw"      -> use signed r as in LD file (no correction)
#   "signed_corrected"-> use signed r corrected for allele swaps after canonicalization
ld_mode <- "abs"

# Threshold rule for lead replication:
#   "abs"    -> require |LD| >= LD_THR
#   "signed" -> require  LD  >= LD_THR  (only positive)
ld_threshold_rule <- "abs"
# --------------------

rep_levels_5 <- c(
  "Exactly replicated",
  "Lead replicated (LD ≥ 0.8)",
  "Lead replicated (LD < 0.8)",
  "Lead replicated (LD missing)",
  "Not replicated"
)

# ==============================================================================
# HELPERS
# ==============================================================================
order_chr <- function(df, chr_col = "CHR") {
  chr_chr <- as.character(df[[chr_col]])
  chr_num <- suppressWarnings(as.numeric(chr_chr))
  ord <- chr_chr[order(chr_num, na.last = TRUE)]
  df[[chr_col]] <- factor(chr_chr, levels = unique(ord))
  df
}

compute_props <- function(df, group_vars) {
  df %>%
    group_by(across(all_of(group_vars))) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup()
}

to_chr_prefix <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- ifelse(is.na(x), NA_character_, x)
  x[x == ""] <- NA_character_
  ifelse(is.na(x), NA_character_,
         ifelse(str_detect(x, "^chr"), x, paste0("chr", x)))
}

extract_pos_from_snpid <- function(x) {
  x <- as.character(x)
  m <- str_match(x, "^(?:chr)?[^:]+:([0-9]+):")
  suppressWarnings(as.numeric(m[,2]))
}

pair_key <- function(a, b, sep = "||") {
  a2 <- as.character(a); b2 <- as.character(b)
  lo <- ifelse(a2 <= b2, a2, b2)
  hi <- ifelse(a2 <= b2, b2, a2)
  paste0(lo, sep, hi)
}

is_valid_snpid <- function(x) {
  x <- as.character(x)
  !is.na(x) & str_detect(x, "^chr[^:]+:[0-9]+:[^:]+:[^:]+$")
}

# ==============================================================================
# Canonical SNPID + flip flag (allele-order invariant ID)
# This function deals with non-harmonized imputed genotype data
# ==============================================================================
canonical_snpid_with_flip <- function(x) {
  x0 <- to_chr_prefix(x)
  m <- stringr::str_match(x0, "^chr([^:]+):([0-9]+):([^:]+):([^:]+)$")
  chr <- m[,2]; pos <- m[,3]; a1 <- m[,4]; a2 <- m[,5]
  
  flip <- ifelse(is.na(a1) | is.na(a2), NA_integer_, as.integer(a1 > a2))
  lo <- ifelse(a1 <= a2, a1, a2)
  hi <- ifelse(a1 <= a2, a2, a1)
  
  id_can <- ifelse(is.na(chr), NA_character_, paste0("chr", chr, ":", pos, ":", lo, ":", hi))
  list(id_can = id_can, flip = flip)
}

# ==============================================================================
# LD MODE FUNCTION (3 options)
# ==============================================================================
apply_ld_mode <- function(r_raw, flip_a, flip_b,
                          mode = c("abs", "signed_raw", "signed_corrected")) {
  mode <- match.arg(mode)
  
  r_raw <- suppressWarnings(as.numeric(r_raw))
  
  if (mode == "abs") return(abs(r_raw))
  if (mode == "signed_raw") return(r_raw)
  
  # signed_corrected
  sgn <- ifelse(is.na(flip_a) | is.na(flip_b), NA_real_,
                ifelse(((flip_a + flip_b) %% 2) == 1, -1, 1))
  r_raw * sgn
}

# ==============================================================================
# LD threshold rule (abs vs signed)
# ==============================================================================
ld_passes_threshold <- function(ld_value, thr = LD_THR, rule = c("abs","signed")) {
  rule <- match.arg(rule)
  ld_value <- as.numeric(ld_value)
  
  out <- rep(FALSE, length(ld_value))
  ok <- !is.na(ld_value)
  
  if (rule == "abs") {
    out[ok] <- abs(ld_value[ok]) >= thr
  } else {
    out[ok] <- ld_value[ok] >= thr
  }
  out
}

# ==============================================================================
# CIS COLLAPSE: best cis per (UNIPROT x REGION)
# ==============================================================================
collapse_cis_by_region <- function(df, region_bp = REGION_BP, anchor_col = ANCHOR_COL) {
  df2 <- df %>%
    mutate(
      CIS_TRANS = ifelse(!is.na(CIS_TRANS) & CIS_TRANS %in% c("cis","Cis","CIS"," cis"), "cis", "trans"),
      CHR_chr    = as.character(CHR),
      anchor_id  = to_chr_prefix(.data[[anchor_col]]),
      anchor_pos = extract_pos_from_snpid(anchor_id),
      region_bin = ifelse(!is.na(anchor_pos), floor(anchor_pos / region_bp), NA_real_),
      region_id  = ifelse(!is.na(CHR_chr) & !is.na(region_bin),
                          paste0("chr", CHR_chr, "_bin", region_bin),
                          NA_character_)
    )
  
  cis_only <- df2 %>% filter(CIS_TRANS == "cis", !is.na(UNIPROT), !is.na(MLOG10P), !is.na(region_id))
  
  cis_best <- cis_only %>%
    group_by(UNIPROT, region_id) %>%
    slice_max(MLOG10P, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  trans_all <- df2 %>% filter(CIS_TRANS == "trans")
  
  bind_rows(cis_best, trans_all) %>%
    select(-CHR_chr, -anchor_id, -anchor_pos)
}

# ==============================================================================
# Needed LD pairs per study (post-collapse)
# ==============================================================================
make_needed_pairs_from_study <- function(path) {
  dat <- read.csv(path, stringsAsFactors = FALSE)
  
  pqtl <- collapse_cis_by_region(dat) %>%
    mutate(CHR = as.character(CHR)) %>%
    filter(CHR %in% as.character(1:22))
  
  # Fallback: if SNPID_EXACT missing/empty -> use SNPID (original study SNP)
  exact_raw_for_ld <- ifelse(is.na(pqtl[[EXACT_ID_COL]]) | pqtl[[EXACT_ID_COL]] == "",
                             pqtl$SNPID,
                             pqtl[[EXACT_ID_COL]])
  
  lead_can  <- canonical_snpid_with_flip(pqtl[[LEAD_ID_COL]])$id_can
  exact_can <- canonical_snpid_with_flip(exact_raw_for_ld)$id_can
  
  pqtl %>%
    mutate(
      lead_id  = lead_can,
      exact_id = exact_can
    ) %>%
    filter(is_valid_snpid(lead_id), is_valid_snpid(exact_id), !is.na(CHR)) %>%
    mutate(ld_key = pair_key(lead_id, exact_id)) %>%
    distinct(CHR, ld_key)
}

# ==============================================================================
# LD loader per chr (filtered) + LD mode
# ==============================================================================
load_ld_chr_filtered <- function(chr, needed_keys, ld_mode = ld_mode) {
  chr <- as.character(chr)
  needed_keys <- unique(needed_keys)
  needed_keys <- needed_keys[!is.na(needed_keys)]
  if (length(needed_keys) == 0) return(tibble(key = character(), LD_r = numeric()))
  
  f <- ld_file_for_chr(chr)
  if (!file.exists(f)) {
    warning("LD file not found for chr ", chr, ": ", f)
    return(tibble(key = character(), LD_r = numeric()))
  }
  
  message("  Loading LD for chr", chr, " (filtered; needed pairs=", length(needed_keys), ")")
  
  ld <- tryCatch(
    fread(cmd = paste("zcat", shQuote(f)), data.table = FALSE, sep = " ", header = TRUE, fill = TRUE, strip.white = TRUE),
    error = function(e) NULL
  )
  if (is.null(ld) || nrow(ld) == 0 || ncol(ld) <= 1) {
    ld <- tryCatch(
      read.table(gzfile(f), header = TRUE, sep = " ", stringsAsFactors = FALSE,
                 quote = "", comment.char = "", fill = TRUE, strip.white = TRUE),
      error = function(e) NULL
    )
  }
  if (is.null(ld) || nrow(ld) == 0 || ncol(ld) <= 1) {
    warning("chr", chr, ": LD file could not be parsed.")
    return(tibble(key = character(), LD_r = numeric()))
  }
  
  names(ld) <- trimws(names(ld))
  if (!all(c("SNP_A","SNP_B","R") %in% names(ld))) {
    stop(sprintf("chr%s: missing SNP_A/SNP_B/R. Found: %s", chr, paste(names(ld), collapse = ", ")))
  }
  
  snpA <- to_chr_prefix(ld$SNP_A)
  snpB <- to_chr_prefix(ld$SNP_B)
  ca <- canonical_snpid_with_flip(snpA)
  cb <- canonical_snpid_with_flip(snpB)
  
  ld2 <- tibble(
    key  = pair_key(ca$id_can, cb$id_can),
    LD_r = apply_ld_mode(ld$R, ca$flip, cb$flip, mode = ld_mode)
  ) %>%
    filter(!is.na(LD_r), !is.na(key)) %>%
    filter(key %in% needed_keys) %>%
    distinct(key, .keep_all = TRUE)
  
  ld2
}

# ==============================================================================
# PreQC presence cache (chr*.bim) + canonicalization
# ==============================================================================
.bim_cache <- new.env(parent = emptyenv())

load_bim_chr_ids_canonical <- function(chr, bfile_dir = bfile_dir) {
  chr <- as.character(chr)
  if (exists(chr, envir = .bim_cache, inherits = FALSE)) {
    return(get(chr, envir = .bim_cache, inherits = FALSE))
  }
  
  bim_path <- file.path(bfile_dir, paste0("chr", chr, ".bim"))
  if (!file.exists(bim_path)) {
    warning("BIM not found for chr ", chr, ": ", bim_path)
    assign(chr, character(0), envir = .bim_cache)
    return(character(0))
  }
  
  bim <- tryCatch(fread(bim_path, header = FALSE, data.table = FALSE, showProgress = FALSE),
                  error = function(e) NULL)
  if (is.null(bim) || ncol(bim) < 2) {
    assign(chr, character(0), envir = .bim_cache)
    return(character(0))
  }
  
  ids_raw <- as.character(bim[[2]])
  ids_can <- canonical_snpid_with_flip(ids_raw)$id_can
  ids_can <- ids_can[!is.na(ids_can)]
  assign(chr, ids_can, envir = .bim_cache)
  ids_can
}

# ==============================================================================
# Add LD + 5-class replication labels
# ==============================================================================
add_ld_and_replication_5class <- function(df, ld_lookup_by_chr, thr_val, ld_thr_val,
                                          ld_threshold_rule = ld_threshold_rule,
                                          bfile_dir = bfile_dir) {
  
  # Use EXACT variant = SNPID (original study SNP) for LD + preQC presence checks
  exact_raw_for_ld <- df[["SNPID"]]
  
  lead_can  <- canonical_snpid_with_flip(df[[LEAD_ID_COL]])$id_can
  exact_can <- canonical_snpid_with_flip(exact_raw_for_ld)$id_can
  
  df2 <- df %>%
    mutate(
      CHR = as.character(CHR),
      lead_id  = lead_can,
      exact_id = exact_can,
      ld_key   = ifelse(is_valid_snpid(lead_id) & is_valid_snpid(exact_id),
                        pair_key(lead_id, exact_id),
                        NA_character_),
      exact_sig = !is.na(MLOG10P_EXACT) & MLOG10P_EXACT >= thr_val,
      lead_sig  = !is.na(MLOG10P_LEAD)  & MLOG10P_LEAD  >= thr_val,
      need_ld   = (!exact_sig) & lead_sig & !is.na(ld_key),
      missing_seqID_in_believe <- is.na(MLOG10P_LEAD) & is.na(MLOG10P_EXACT) 
    )
  
  # ---- Join LD only for rows that need it ----
  df2 <- df2 %>%
    group_by(CHR) %>%
    group_modify(function(.x, .g) {
      chr <- as.character(.g$CHR[[1]])
      ld_tbl <- ld_lookup_by_chr[[chr]]
      
      LD_r_out <- rep(NA_real_, nrow(.x))
      if (is.null(ld_tbl) || nrow(ld_tbl) == 0) {
        .x$LD_r <- LD_r_out
        return(.x)
      }
      
      idx <- which(.x$need_ld)
      if (length(idx) == 0) {
        .x$LD_r <- LD_r_out
        return(.x)
      }
      
      tmp <- .x[idx, , drop = FALSE] %>%
        left_join(ld_tbl, by = c("ld_key" = "key"))  # ld_tbl has LD_r
      
      LD_r_out[idx] <- tmp$LD_r
      .x$LD_r <- LD_r_out
      .x
    }) %>%
    ungroup()
  
  # ---- preQC presence only where need_ld is TRUE ----
  df2 <- df2 %>%
    group_by(CHR) %>%
    group_modify(function(.x, .g) {
      chr <- as.character(.g$CHR[[1]])
      ids <- load_bim_chr_ids_canonical(chr, bfile_dir = bfile_dir)
      
      .x$lead_in_preQC  <- NA
      .x$exact_in_preQC <- NA
      if (length(ids) == 0) return(.x)
      
      idx <- which(.x$need_ld)
      if (length(idx) == 0) return(.x)
      
      .x$lead_in_preQC[idx]  <- .x$lead_id[idx]  %in% ids
      .x$exact_in_preQC[idx] <- .x$exact_id[idx] %in% ids
      .x
    }) %>%
    ungroup()
  
  df2 <- df2 %>%
    mutate(
      ld_missing_due_to_preQC = case_when(
        !need_ld ~ NA,  # we don't care outside lead-candidates
        is.na(lead_in_preQC) | is.na(exact_in_preQC) ~ NA,
        !(lead_in_preQC & exact_in_preQC) ~ TRUE,
        TRUE ~ FALSE
      ),
      
      ld_bucket = case_when(
        !is.na(LD_r) & ld_passes_threshold(LD_r, thr = ld_thr_val, rule = ld_threshold_rule) ~ "LD_ge",
        !is.na(LD_r) & !ld_passes_threshold(LD_r, thr = ld_thr_val, rule = ld_threshold_rule) ~ "LD_lt",
        is.na(LD_r) & ld_missing_due_to_preQC == TRUE ~ "LD_missing",
        TRUE ~ "LD_lt"  # LD not in file but variants present -> effectively < threshold
      ),
      
      replicated_5class = case_when(
        exact_sig ~ "Exactly replicated",
        !exact_sig & lead_sig & ld_bucket == "LD_ge"      ~ "Lead replicated (LD ≥ 0.6)",
        !exact_sig & lead_sig & ld_bucket == "LD_lt"      ~ "Lead replicated (LD < 0.6)",
        !exact_sig & lead_sig & ld_bucket == "LD_missing" ~ "Lead replicated (LD missing)",
        TRUE ~ "Not replicated"
      ),
      replicated_5class = factor(replicated_5class, levels = rep_levels_5)
    )
  
  df2
}
# ==============================================================================
# Per-study processor
# ==============================================================================
process_study <- function(study, path, ld_lookup_by_chr, thr_val = thr, ld_thr_val = LD_THR,
                          ld_threshold_rule = ld_threshold_rule,
                          bfile_dir = bfile_dir) {
  message("Processing: ", study)
  dat <- read.csv(path, stringsAsFactors = FALSE)
  
  needed <- c("UNIPROT","CIS_TRANS","MLOG10P","MLOG10P_LEAD","MLOG10P_EXACT","CHR",
              "SNPID", LEAD_ID_COL, EXACT_ID_COL)
  miss <- setdiff(needed, names(dat))
  if (length(miss) > 0) stop(sprintf("Study '%s' missing columns: %s", study, paste(miss, collapse=", ")))

  pqtl <- collapse_cis_by_region(dat) %>%
    mutate(CHR = as.character(CHR)) %>%
    filter(CHR %in% as.character(1:22))
  pqtl <- add_ld_and_replication_5class(pqtl, ld_lookup_by_chr, thr_val, ld_thr_val,
                                        ld_threshold_rule = ld_threshold_rule,
                                        bfile_dir = bfile_dir)
  pqtl <- pqtl %>% mutate(CHR = as.character(CHR))
  
  rep_prop_ct <- pqtl %>%
    count(CIS_TRANS, replicated_5class, name="n") %>%
    complete(CIS_TRANS, replicated_5class, fill=list(n=0)) %>%
    compute_props(group_vars = c("CIS_TRANS")) %>%
    mutate(study = study, .before=1)
  
  rep_prop_overall <- pqtl %>%
    count(replicated_5class, name="n") %>%
    complete(replicated_5class, fill=list(n=0)) %>%
    mutate(CIS_TRANS="Overall") %>%
    compute_props(group_vars = c("CIS_TRANS")) %>%
    mutate(study = study, .before=1)
  
  list(study=study, pqtl=pqtl, rep_prop_ct=rep_prop_ct, rep_prop_overall=rep_prop_overall)
}

# ==============================================================================
# RUN ALL STUDIES (forward literature review: literature study signals -> BELIEVE replication)
# ==============================================================================
study_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
study_dirs <- study_dirs[grepl("^pqtl_", basename(study_dirs))]

studies <- tibble(
  study = basename(study_dirs),
  path  = file.path(study_dirs, paste0(basename(study_dirs), "_hdsc_believe.csv"))
) %>%
  mutate(study_clean = str_remove(study, "^pqtl_"))

# edit if needed
studies_filt <- studies %>%
 filter(study_clean %in% c("sun","sun_ukb", "sun_ukb_csa", "interval_chris_meta","fenland","QBB","JHS"))
# c("kora", "ElderlyEU","AGES-65","aric_EA","Brain","CSF","QMDiab","IBD_Europe","Plasma"))
# studies_filt <- studies %>%
#   filter(study_clean %in% c("JHS"))

# PASS 1: needed LD keys
message("Collecting needed LD pairs across studies (post-collapse)...")
needed_pairs_all <- map_dfr(studies_filt$path, make_needed_pairs_from_study) %>%
  distinct(CHR, ld_key)
needed_keys_by_chr <- split(needed_pairs_all$ld_key, needed_pairs_all$CHR)

# PASS 2: LD lookup
message("Building per-chromosome LD lookup tables (filtered)...")
ld_lookup_by_chr <- list()
for (chr in names(needed_keys_by_chr)) {
  ld_lookup_by_chr[[chr]] <- load_ld_chr_filtered(chr, needed_keys_by_chr[[chr]], ld_mode = ld_mode)
}

# PASS 3: process studies
results_list <- pmap(
  studies_filt,
  ~ process_study(study = ..1, path = ..2,
                  ld_lookup_by_chr = ld_lookup_by_chr,
                  thr_val = thr, ld_thr_val = LD_THR,
                  ld_threshold_rule = ld_threshold_rule,
                  bfile_dir = bfile_dir)
)
names(results_list) <- studies_filt$study

# Combined tables for plots
all_rep_prop_ct      <- bind_rows(map(results_list, "rep_prop_ct"))
all_rep_prop_overall <- bind_rows(map(results_list, "rep_prop_overall"))

# ==============================================================================
# COUNTS YOU ASKED FOR:
# how many lead-replicated have LD < 0.6 vs missing due to preQC absence
# ==============================================================================
lead_ld_breakdown <- map_dfr(names(results_list), function(st) {
  df <- results_list[[st]]$pqtl
  
  lead_sig_exact_not <- df %>%
    filter(!is.na(MLOG10P_LEAD) & MLOG10P_LEAD >= thr,
           is.na(MLOG10P_EXACT) | MLOG10P_EXACT < thr)
  
  tibble(
    study = st,
    n_lead_sig_exact_not = nrow(lead_sig_exact_not),
    n_lead_rep_ld_ge_08  = sum(lead_sig_exact_not$replicated_5class == "Lead replicated (LD ≥ 0.6)", na.rm = TRUE),
    n_lead_rep_ld_lt_08  = sum(lead_sig_exact_not$replicated_5class == "Lead replicated (LD < 0.6)", na.rm = TRUE),
    n_lead_rep_ld_missing_preQC = sum(lead_sig_exact_not$replicated_5class == "Lead replicated (LD missing)", na.rm = TRUE)
  )
}) %>% arrange(desc(n_lead_sig_exact_not))

print(lead_ld_breakdown)

# ==============================================================================
# PLOTS (5 classes)
# ==============================================================================
cb_pal <- c(
  "Exactly replicated" = "#1B7837",  
  "Lead replicated (LD ≥ 0.6)" = "green",  #
  "Lead replicated (LD < 0.6)" = "yellow",  
  "Lead replicated (LD missing)" = "#E69F00", 
  "Not replicated" = "#2166AC"  # black
)

p_rep_overall <- ggplot(all_rep_prop_overall,
                        aes(x = study, y = prop, fill = replicated_5class)) +
  geom_col() +
  scale_fill_manual(values = cb_pal) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = NULL, y = "Percentage", fill = "Replication",
       title = "Forward replication (Overall)") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        legend.position = "bottom")
print(p_rep_overall)

p_rep_ct <- ggplot(all_rep_prop_ct,
                   aes(x = CIS_TRANS, y = prop, fill = replicated_5class)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = cb_pal) +
  scale_y_continuous(labels = scales::percent_format()) +
  facet_wrap(~ study) +
  labs(x = NULL, y = "Percentage", fill = "Replication",
       title = "Forward replication (cis vs trans)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")
print(p_rep_ct)

# ==============================================================================
# OPTIONAL: sanity check on LD NA reasons for one study of the list
# ==============================================================================
st <- names(results_list)[1]
df <- results_list[[st]]$pqtl
print(table(is.na(df$LD_r), df$ld_missing_due_to_preQC, useNA = "ifany"))
