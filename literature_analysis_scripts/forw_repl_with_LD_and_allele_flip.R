# ==============================================================================
# FORWARD LITERATURE REVIEW PIPELINE 
# - 8 replication classes
# - LD mode: abs / signed_raw / signed_corrected
# - LD threshold rule: abs (|LD|>=thr) / signed (LD>=thr)
# - Distinguishes LD missing due to variants absent in preQC (.bim)
#   vs LD absent from LD file (filtered) while variants exist in preQC
# - Expanded replication classes (concordant/discordant)
# - Lead replication uses |r| threshold, concordance uses sign(r)
# - Adds "Missing protein" when UNIPROT_MATCH=="" and seqID missing in BELIEVE
# - LD threshold is set ONCE at the top (LD_R_THR)
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
# SETTINGS (EDIT HERE)
# ==============================================================================
thr        <- -log10(5e-8)

LD_R_THR   <- 0.8       # threshold on r (NOT r^2). e.g. 0.8 -> r^2 >= 0.64
REGION_BP  <- 1e6
ANCHOR_COL <- "SNPID_LEAD"

base_dir   <- "/exchange/healthds/pQTL/BELIEVE/gwasstudio_output"
bfile_dir  <- "/scratch/laura.bondi/BELIEVE_bfile"

ld_base_dir <- "/scratch/laura.bondi/PLINK_LD_RESUME_FINAL"
ld_file_for_chr <- function(chr) file.path(ld_base_dir, paste0("chr", chr), paste0("chr", chr, "_r.ld.gz"))

LEAD_ID_COL  <- "SNPID_LEAD"  # BELIEVE lead
EXACT_ID_COL <- "SNPID"       # published exact SNP (always)

# LD sign handling for concordance (recommended: signed_corrected)
ld_mode_signed <- "signed_corrected"  # "signed_raw" or "signed_corrected"

# ==============================================================================
# CLASS LEVELS (expanded)
# ==============================================================================
rep_levels <- c(
  "Exactly replicated (Concordant)",
  "Exactly replicated (Discordant)",
  sprintf("Lead replicated (|LD| ≥ %.2f; Concordant)", LD_R_THR),
  sprintf("Lead replicated (|LD| ≥ %.2f; Discordant)", LD_R_THR),
  sprintf("Lead replicated (|LD| < %.2f)", LD_R_THR),
  "Lead replicated (LD missing)",
  "Missing protein",
  "Not replicated"
)

# ==============================================================================
# HELPERS
# ==============================================================================
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

apply_ld_signed <- function(r_raw, flip_a, flip_b, mode = c("signed_raw","signed_corrected")) {
  mode <- match.arg(mode)
  r_raw <- suppressWarnings(as.numeric(r_raw))
  if (mode == "signed_raw") return(r_raw)
  
  sgn <- ifelse(is.na(flip_a) | is.na(flip_b), NA_real_,
                ifelse(((flip_a + flip_b) %% 2) == 1, -1, 1))
  r_raw * sgn
}

beta_sign <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.na(x) | x == 0, NA_integer_, sign(x))
}

exact_concordant_flag <- function(beta_exact, beta_pub) {
  s_exact <- beta_sign(beta_exact)
  s_pub   <- beta_sign(beta_pub)
  ifelse(is.na(s_exact) | is.na(s_pub), NA, (s_exact == s_pub))
}

lead_concordant_flag <- function(beta_lead, beta_pub, ld_r_signed) {
  s_lead <- beta_sign(beta_lead)
  s_pub  <- beta_sign(beta_pub)
  s_ld   <- beta_sign(ld_r_signed)
  ifelse(is.na(s_lead) | is.na(s_pub) | is.na(s_ld), NA, (s_lead == s_pub * s_ld))
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
# PASS 1: Needed LD keys (only where LD is needed)
# ==============================================================================
make_needed_pairs_from_study <- function(path) {
  dat <- read.csv(path, stringsAsFactors = FALSE)
  
  pqtl <- collapse_cis_by_region(dat) %>%
    mutate(CHR = as.character(CHR)) %>%
    filter(CHR %in% as.character(1:22))
  
  lead_can  <- canonical_snpid_with_flip(pqtl[[LEAD_ID_COL]])$id_can
  exact_can <- canonical_snpid_with_flip(pqtl[[EXACT_ID_COL]])$id_can
  
  pqtl %>%
    mutate(
      lead_id = lead_can,
      exact_id = exact_can,
      exact_sig = !is.na(MLOG10P_EXACT) & MLOG10P_EXACT >= thr,
      lead_sig  = !is.na(MLOG10P_LEAD)  & MLOG10P_LEAD  >= thr
    ) %>%
    filter(!exact_sig, lead_sig) %>%  # LD only needed here
    filter(is_valid_snpid(lead_id), is_valid_snpid(exact_id), !is.na(CHR)) %>%
    mutate(ld_key = pair_key(lead_id, exact_id)) %>%
    distinct(CHR, ld_key)
}

# ==============================================================================
# PASS 2: Load LD for one chr (filtered) -> returns signed+abs r
# ==============================================================================
load_ld_chr_filtered <- function(chr, needed_keys, ld_mode_signed = ld_mode_signed) {
  chr <- as.character(chr)
  needed_keys <- unique(needed_keys)
  needed_keys <- needed_keys[!is.na(needed_keys)]
  if (length(needed_keys) == 0) return(tibble(key = character(), LD_r_signed = numeric(), LD_r_abs = numeric()))
  
  f <- ld_file_for_chr(chr)
  if (!file.exists(f)) {
    warning("LD file not found for chr ", chr, ": ", f)
    return(tibble(key = character(), LD_r_signed = numeric(), LD_r_abs = numeric()))
  }
  
  message("  Loading LD for chr", chr, " (filtered; needed pairs=", length(needed_keys), ")")
  
  ld <- tryCatch(
    fread(cmd = paste("zcat", shQuote(f)),
          data.table = FALSE, sep = " ", header = TRUE, fill = TRUE, strip.white = TRUE),
    error = function(e) NULL
  )
  if (is.null(ld) || nrow(ld) == 0 || ncol(ld) <= 1) {
    ld <- tryCatch(
      read.table(gzfile(f), header = TRUE, sep = "", stringsAsFactors = FALSE,
                 quote = "", comment.char = "", fill = TRUE, strip.white = TRUE),
      error = function(e) NULL
    )
  }
  if (is.null(ld) || nrow(ld) == 0 || ncol(ld) <= 1) {
    warning("chr", chr, ": LD file could not be parsed.")
    return(tibble(key = character(), LD_r_signed = numeric(), LD_r_abs = numeric()))
  }
  
  names(ld) <- trimws(names(ld))
  if (!all(c("SNP_A","SNP_B","R") %in% names(ld))) {
    stop(sprintf("chr%s: missing SNP_A/SNP_B/R. Found: %s", chr, paste(names(ld), collapse = ", ")))
  }
  
  snpA <- to_chr_prefix(ld$SNP_A)
  snpB <- to_chr_prefix(ld$SNP_B)
  ca <- canonical_snpid_with_flip(snpA)
  cb <- canonical_snpid_with_flip(snpB)
  
  r_signed <- apply_ld_signed(ld$R, ca$flip, cb$flip, mode = ld_mode_signed)
  
  tibble(
    key = pair_key(ca$id_can, cb$id_can),
    LD_r_signed = suppressWarnings(as.numeric(r_signed)),
    LD_r_abs    = abs(suppressWarnings(as.numeric(r_signed)))
  ) %>%
    filter(!is.na(key), !is.na(LD_r_signed)) %>%
    filter(key %in% needed_keys) %>%
    distinct(key, .keep_all = TRUE)
}

# ==============================================================================
# PreQC presence cache (chr*.bim)
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
  
  ids_raw <- trimws(as.character(bim[[2]]))
  ids_can <- canonical_snpid_with_flip(ids_raw)$id_can
  ids_can <- ids_can[!is.na(ids_can)]
  assign(chr, ids_can, envir = .bim_cache)
  ids_can
}

# ==============================================================================
# PASS 3: Add LD + expanded classes + concordance flags
# ==============================================================================
add_ld_and_replication_expanded <- function(df, ld_lookup_by_chr, thr_val, ld_thr_val = LD_R_THR,
                                            bfile_dir = bfile_dir) {
  
  lead_can  <- canonical_snpid_with_flip(df[[LEAD_ID_COL]])$id_can
  exact_can <- canonical_snpid_with_flip(df[[EXACT_ID_COL]])$id_can
  
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
      
      # seqID not found in BELIEVE -> no p-values at all
      missing_seqID_in_believe = is.na(MLOG10P_LEAD) & is.na(MLOG10P_EXACT),
      
      # Missing protein definition you requested
      missing_protein = missing_seqID_in_believe & (trimws(UNIPROT_MATCH) == "")
    )
  
  # Join LD only where need_ld is TRUE
  df2 <- df2 %>%
    group_by(CHR) %>%
    group_modify(function(.x, .g) {
      chr <- as.character(.g$CHR[[1]])
      ld_tbl <- ld_lookup_by_chr[[chr]]
      
      LD_signed <- rep(NA_real_, nrow(.x))
      LD_abs    <- rep(NA_real_, nrow(.x))
      
      if (!is.null(ld_tbl) && nrow(ld_tbl) > 0) {
        idx <- which(.x$need_ld)
        if (length(idx) > 0) {
          tmp <- .x[idx, , drop = FALSE] %>% left_join(ld_tbl, by = c("ld_key" = "key"))
          LD_signed[idx] <- tmp$LD_r_signed
          LD_abs[idx]    <- tmp$LD_r_abs
        }
      }
      
      .x$LD_r_signed <- LD_signed
      .x$LD_r_abs    <- LD_abs
      .x
    }) %>%
    ungroup()
  
  # preQC presence only where need_ld is TRUE
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
      # only meaningful where need_ld is TRUE
      ld_missing_due_to_preQC = case_when(
        !need_ld ~ NA,
        is.na(lead_in_preQC) | is.na(exact_in_preQC) ~ NA,
        !(lead_in_preQC & exact_in_preQC) ~ TRUE,
        TRUE ~ FALSE
      ),
      
      # LD bucket uses ABS(r)
      ld_bucket = case_when(
        !is.na(LD_r_abs) & (LD_r_abs >= ld_thr_val) ~ "LD_ge",
        !is.na(LD_r_abs) & (LD_r_abs <  ld_thr_val) ~ "LD_lt",
        is.na(LD_r_abs) & ld_missing_due_to_preQC == TRUE ~ "LD_missing",
        TRUE ~ "LD_lt"   # not in LD file but variants exist -> effectively below threshold
      ),
      
      # concordance flags
      exact_concordant = ifelse(exact_sig, exact_concordant_flag(BETA_EXACT, BETA), NA),
      
      lead_concordant  = ifelse((!exact_sig) & lead_sig & (ld_bucket == "LD_ge"),
                                lead_concordant_flag(BETA_LEAD, BETA, LD_r_signed), NA),
      
      replicated_class = case_when(
        # New class first
        missing_protein ~ "Missing protein",
        
        # Exact replicated split
        exact_sig & exact_concordant == TRUE  ~ "Exactly replicated (Concordant)",
        exact_sig & exact_concordant == FALSE ~ "Exactly replicated (Discordant)",
        exact_sig & is.na(exact_concordant)   ~ "Exactly replicated (Concordant)",
        
        # Lead replicated with high |LD| split
        (!exact_sig) & lead_sig & (ld_bucket == "LD_ge") & lead_concordant == TRUE  ~ sprintf("Lead replicated (|LD| ≥ %.2f; Concordant)", ld_thr_val),
        (!exact_sig) & lead_sig & (ld_bucket == "LD_ge") & lead_concordant == FALSE ~ sprintf("Lead replicated (|LD| ≥ %.2f; Discordant)", ld_thr_val),
        (!exact_sig) & lead_sig & (ld_bucket == "LD_ge") & is.na(lead_concordant)   ~ sprintf("Lead replicated (|LD| ≥ %.2f; Concordant)", ld_thr_val),
        
        # Other lead categories
        (!exact_sig) & lead_sig & (ld_bucket == "LD_lt")      ~ sprintf("Lead replicated (|LD| < %.2f)", ld_thr_val),
        (!exact_sig) & lead_sig & (ld_bucket == "LD_missing") ~ "Lead replicated (LD missing)",
        
        TRUE ~ "Not replicated"
      ),
      replicated_class = factor(replicated_class, levels = rep_levels)
    )
  
  df2
}

# ==============================================================================
# Per-study processor
# ==============================================================================
process_study <- function(study, path, ld_lookup_by_chr, thr_val = thr, ld_thr_val = LD_R_THR,
                          bfile_dir = bfile_dir) {
  
  message("Processing: ", study)
  dat <- read.csv(path, stringsAsFactors = FALSE)
  
  needed <- c("UNIPROT","UNIPROT_MATCH","CIS_TRANS","MLOG10P","MLOG10P_LEAD","MLOG10P_EXACT","CHR",
              "SNPID","BETA","BETA_LEAD","BETA_EXACT", LEAD_ID_COL, EXACT_ID_COL)
  miss <- setdiff(needed, names(dat))
  if (length(miss) > 0) stop(sprintf("Study '%s' missing columns: %s", study, paste(miss, collapse=", ")))
  
  pqtl <- collapse_cis_by_region(dat) %>%
    mutate(CHR = as.character(CHR)) %>%
    filter(CHR %in% as.character(1:22))
  
  pqtl <- add_ld_and_replication_expanded(
    pqtl,
    ld_lookup_by_chr = ld_lookup_by_chr,
    thr_val = thr_val,
    ld_thr_val = ld_thr_val,
    bfile_dir = bfile_dir
  )
  
  rep_prop_ct <- pqtl %>%
    count(CIS_TRANS, replicated_class, name="n") %>%
    complete(CIS_TRANS, replicated_class, fill=list(n=0)) %>%
    compute_props(group_vars = c("CIS_TRANS")) %>%
    mutate(study = study, .before=1)
  
  rep_prop_overall <- pqtl %>%
    count(replicated_class, name="n") %>%
    complete(replicated_class, fill=list(n=0)) %>%
    mutate(CIS_TRANS="Overall") %>%
    compute_props(group_vars = c("CIS_TRANS")) %>%
    mutate(study = study, .before=1)
  
  list(study=study, pqtl=pqtl, rep_prop_ct=rep_prop_ct, rep_prop_overall=rep_prop_overall)
}

# ==============================================================================
# RUN STUDIES
# ==============================================================================
study_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
study_dirs <- study_dirs[grepl("^pqtl_", basename(study_dirs))]

studies <- tibble(
  study = basename(study_dirs),
  path  = file.path(study_dirs, paste0(basename(study_dirs), "_hdsc_believe.csv"))
) %>%
  mutate(study_clean = str_remove(study, "^pqtl_"))

# Edit this selection
studies_filt <- studies %>%
  filter(study_clean %in% c("interval_chris_meta","QBB"))

# PASS 1: needed LD keys
message("Collecting needed LD pairs across studies (post-collapse)...")
needed_pairs_all <- map_dfr(studies_filt$path, make_needed_pairs_from_study) %>%
  distinct(CHR, ld_key)
needed_keys_by_chr <- split(needed_pairs_all$ld_key, needed_pairs_all$CHR)

# PASS 2: LD lookup
message("Building per-chromosome LD lookup tables (filtered)...")
ld_lookup_by_chr <- list()
for (chr in names(needed_keys_by_chr)) {
  ld_lookup_by_chr[[chr]] <- load_ld_chr_filtered(chr, needed_keys_by_chr[[chr]], ld_mode_signed = ld_mode_signed)
}

# PASS 3: process studies
results_list <- pmap(
  studies_filt,
  ~ process_study(study = ..1, path = ..2,
                  ld_lookup_by_chr = ld_lookup_by_chr,
                  thr_val = thr, ld_thr_val = LD_R_THR,
                  bfile_dir = bfile_dir)
)
names(results_list) <- studies_filt$study

all_rep_prop_ct      <- bind_rows(map(results_list, "rep_prop_ct"))
all_rep_prop_overall <- bind_rows(map(results_list, "rep_prop_overall"))

# ==============================================================================
# COUNTS: lead replicated break-down (includes concordance split)
# ==============================================================================
lead_ld_breakdown <- map_dfr(names(results_list), function(st) {
  df <- results_list[[st]]$pqtl
  
  lead_candidates <- df %>%
    filter(!is.na(MLOG10P_LEAD) & MLOG10P_LEAD >= thr,
           is.na(MLOG10P_EXACT) | MLOG10P_EXACT < thr)
  
  tibble(
    study = st,
    n_lead_candidates = nrow(lead_candidates),
    n_lead_ld_ge_conc = sum(lead_candidates$replicated_class == sprintf("Lead replicated (|LD| ≥ %.2f; Concordant)", LD_R_THR), na.rm = TRUE),
    n_lead_ld_ge_disc = sum(lead_candidates$replicated_class == sprintf("Lead replicated (|LD| ≥ %.2f; Discordant)", LD_R_THR), na.rm = TRUE),
    n_lead_ld_lt      = sum(lead_candidates$replicated_class == sprintf("Lead replicated (|LD| < %.2f)", LD_R_THR), na.rm = TRUE),
    n_lead_ld_missing = sum(lead_candidates$replicated_class == "Lead replicated (LD missing)", na.rm = TRUE)
  )
})

print(lead_ld_breakdown)

# ==============================================================================
# PLOTS (expanded classes)
# ==============================================================================
cb_pal <- c(
  "Exactly replicated (Concordant)" = "#009E73",
  "Exactly replicated (Discordant)" = "#D55E00",
  sprintf("Lead replicated (|LD| ≥ %.2f; Concordant)", LD_R_THR) = "#0072B2",
  sprintf("Lead replicated (|LD| ≥ %.2f; Discordant)", LD_R_THR) = "#CC79A7",
  sprintf("Lead replicated (|LD| < %.2f)", LD_R_THR) = "#E69F00",
  "Lead replicated (LD missing)" = "#999999",
  "Missing protein" = "#000000",
  "Not replicated" = "#56B4E9"
)

p_rep_overall <- ggplot(all_rep_prop_overall,
                        aes(x = study, y = prop, fill = replicated_class)) +
  geom_col() +
  scale_fill_manual(values = cb_pal, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = NULL, y = "Percentage", fill = "Replication",
       title = "Forward replication (Overall)") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        legend.position = "bottom")
print(p_rep_overall)

p_rep_ct <- ggplot(all_rep_prop_ct,
                   aes(x = CIS_TRANS, y = prop, fill = replicated_class)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = cb_pal, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format()) +
  facet_wrap(~ study) +
  labs(x = NULL, y = "Percentage", fill = "Replication",
       title = "Forward replication (cis vs trans)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")
print(p_rep_ct)
