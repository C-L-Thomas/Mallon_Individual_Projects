# ============================================================
# Bumblebee methylation structure and reproductive status
# ------------------------------------------------------------
# Molecular Ecology revision - k=10 CENTROID-DISTANCE, OMNIBUS script
#
# ------------------------------------------------------------
# "_omnibus" VARIANT NOTE
# Built from bombus_univariate_distance_reanalysis_centroid_k10_2.R. Changes:
#
#   1. REMOVED the k=2 / k=3 comparison (PART B1's fit_other_k() and
#      compare_k_all). That question is settled; this script is k=10 only.
#   2. STANDALONE for all three colony-label variants. The old script only
#      ran the "corrected" labels. This one loops over corrected / excluded
#      / uncorrected, each refit from raw data (not sliced from a shared
#      fit), and assembles the three-way comparison table itself
#      (manuscript Supplementary Table 15) - no need to also run the sibling
#      _samples_removed_ / _uncorrected_ scripts to get that table.
#   3. RE-ADDED the old script's PART 3a (eigenvalue/variance-explained
#      table + PC1 permutation nulls, classical & robust), which had been
#      dropped entirely from the k10_2 draft. Gated behind RUN_PART3A at
#      the top, since it's slow (refits PCA up to ~50 times on the full
#      ~285k-CpG matrix). Runs on the primary (corrected) variant only.
#   4. RE-ADDED the UMAP colony-clustering figure (main + common-axis
#      supplementary version), also dropped from the k10_2 draft. Now
#      naturally uses all 10 PCs from the k=10 fit (the original script's
#      UMAP was capped at 5 PCs, a side-effect of that script's k_extract=5
#      PCA - not a deliberate choice).
#   5. FIXED the predicted-probability curve (Figure 4/5 panel A). The
#      k10_2 draft plotted dist_k_ref directly (an earlier, superseded
#      "absolute distance + arrows" style carried over from mid-development
#      of the original script). The manuscript describes, and the original
#      script's *final* combined figure actually used, distance expressed
#      *relative to the matched pair partner within the same stratum*
#      (rel_dist_ref <- PC1_dist_ref - rev(PC1_dist_ref), which works
#      because each stratum has exactly one sterile + one reproductive
#      individual). That paired version is what's built here
#      (add_relative_distance() / rel_dist_k_ref), generalised to the k=10
#      centroid distance. The absolute/arrow version is not reproduced.
#   6. ADDED the PC1+PC2 linear-model check the manuscript cites ("A model
#      including the first two robust methylation principal components as
#      linear predictors provided no evidence of association...
#      chisq2=0.74, p=0.70"). That number was never in either centroid
#      script and, in the old univariate script, came from a k_extract=5
#      PCA fit. It is recomputed here fresh from THIS script's own k=10 PCA
#      fit for every label variant - so do not be surprised if the primary
#      (corrected) variant's number differs from 0.74/0.70 in the current
#      manuscript draft; that old number is not valid for k=10 and should
#      be replaced with whatever this script reports.
#   7. EXTENDED PART 5's robust-PCA permutation null from PC1-only to all
#      k_dist=10 components: a proper parallel-analysis-style test, where
#      each permutation replicate is fit at k=10 (one PcaHubert call per
#      replicate, same as before - NOT one call per component), so every
#      component's observed variance share is compared against its own
#      null distribution rather than just PC1's. The separate PC1-only
#      robust test is gone - PC1's result is just the first row of the new
#      per-component table. The classical (non-robust) PCA permutation
#      null is now OFF by default (RUN_CLASSICAL_PERM_NULL toggle): it was
#      a robust-vs-classical cross-check on the spectrum shape, not part of
#      the core "why k=10" argument, and PART 8 already cross-checks
#      classical vs. robust PCA on the model result itself. The code is
#      kept, just gated off, in case you want it back for the SI.
#
# All outputs from this script carry an "_omnibus" suffix.
# ------------------------------------------------------------

# ============================================================
# PART 0: LIBRARIES, TOGGLES & SETTINGS
# ============================================================

library(data.table)
library(survival)
library(splines)
library(ggplot2)
library(dplyr)
library(rrcov)
library(uwot)
library(patchwork)

## ---- Toggles ----
RUN_PART3A      <- TRUE   # eigenvalue table + per-component robust permutation nulls
                           # (SLOW: refits robust PCA ~49 times on the full ~285k-CpG
                           # matrix). Set FALSE to skip straight to the model results.
RUN_CLASSICAL_PERM_NULL <- FALSE   # optional cross-check nested under RUN_PART3A: a
                           # classical (non-robust) PCA, PC1-only permutation null
                           # (199 reps of plain prcomp - fast, not the bottleneck).
                           # Off by default; see header note, point 7, for why. Set
                           # TRUE to bring it back.
RUN_SENSITIVITY <- TRUE   # 5-variant sensitivity sweep (transform / PCA method /
                           # imputation), primary variant only (SLOW: 5 more PCA
                           # refits). Set FALSE to skip. Added for the same reason as
                           # RUN_PART3A even though only RUN_PART3A was requested
                           # explicitly - remove this toggle if you'd rather it always run.

LABEL_VARIANTS  <- c("corrected", "excluded", "uncorrected")
PRIMARY_VARIANT <- "corrected"   # drives PART 5 (PART 3a), PART 6 (UMAP), and PART 9 (figures)

k_dist <- 10   # fixed for this script - see header note. PCA is always refit fresh
               # at exactly this k, never sliced from a larger existing fit.

DISPUTED_SAMPLES <- c("N6B3-3", "N6B3-4", "N7B2-1", "N7B2-2")

out_dir <- "/Users/emb3/Library/CloudStorage/OneDrive-UniversityofLeicester/projects/Eamonn_Bumblebee/copy_for_paper/revision/outputs"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# PART 1: RAW DATA IMPORT (unmodified - label variants are applied downstream,
# per-variant, in PART 3)
# ============================================================

C_raw <- read.table(
  "/Users/emb3/Library/CloudStorage/OneDrive-UniversityofLeicester/projects/Eamonn_Bumblebee/methylated_counts_matrix.tsv",
  header = TRUE, row.names = 1, sep = "\t", check.names = FALSE
)

N_raw <- read.table(
  "/Users/emb3/Library/CloudStorage/OneDrive-UniversityofLeicester/projects/Eamonn_Bumblebee/coverage_counts_matrix.tsv",
  header = TRUE, row.names = 1, sep = "\t", check.names = FALSE
)

covar_raw <- fread(
  "/Users/emb3/Library/CloudStorage/OneDrive-UniversityofLeicester/projects/Eamonn_Bumblebee/covariates.tsv"
)

if (!("treatment" %in% names(covar_raw))) {
  if ("group" %in% names(covar_raw)) {
    covar_raw[, treatment := group]
  } else {
    stop("Neither 'treatment' nor 'group' found in covariates.tsv. Columns are: ",
         paste(names(covar_raw), collapse = ", "))
  }
}

# As loaded, covar_raw$colony is the AS-RECORDED (uncorrected) state: sample-name
# prefix matches colony (N6B3-3/N6B3-4 -> "N6", N7B2-1/N7B2-2 -> "N7"). The
# genetic-evidence correction (see manuscript Methods) swaps these four to the
# opposite colony; "excluded" drops them entirely rather than picking a label.

# ============================================================
# PART 2: HELPER FUNCTIONS
# ============================================================

## ---- 2.1: build the covariate table for one label variant ----
build_covar_variant <- function(covar_raw, variant) {
  cv <- data.table::copy(covar_raw)
  if (variant == "corrected") {
    cv[sample %in% c("N6B3-3", "N6B3-4"), colony := "N7"]
    cv[sample %in% c("N7B2-1", "N7B2-2"), colony := "N6"]
  } else if (variant == "uncorrected") {
    # as-recorded: no change
  } else if (variant == "excluded") {
    cv <- cv[!(sample %in% DISPUTED_SAMPLES)]
  } else {
    stop("Unknown label variant: ", variant)
  }
  cv
}

## ---- 2.2: align C/N to a covariate table; build methylation-proportion matrices ----
build_meth_matrices <- function(C_raw, N_raw, cv) {
  common <- intersect(colnames(C_raw), cv$sample)
  C <- C_raw[, common, drop = FALSE]
  N <- N_raw[, common, drop = FALSE]
  cv <- cv[match(common, cv$sample)]
  stopifnot(identical(colnames(C), cv$sample))
  stopifnot(identical(colnames(N), cv$sample))

  C_mat <- as.matrix(C); storage.mode(C_mat) <- "double"
  N_mat <- as.matrix(N); storage.mode(N_mat) <- "double"

  meth_prop <- C_mat / N_mat
  meth_prop[!is.finite(meth_prop)] <- NA
  meth_prop_raw <- meth_prop   # pre-imputation snapshot, needed for complete-case sensitivity

  if (anyNA(meth_prop)) {
    rm <- rowMeans(meth_prop, na.rm = TRUE)
    na_rows <- which(rowSums(is.na(meth_prop)) > 0)
    for (i in na_rows) meth_prop[i, is.na(meth_prop[i, ])] <- rm[i]
  }

  list(cv = cv, C_mat = C_mat, N_mat = N_mat, meth_prop = meth_prop, meth_prop_raw = meth_prop_raw)
}

## ---- 2.3: fit robust PCA at exactly k = k_dist and build `dat` ----
fit_k10_and_build_dat <- function(meth_prop, cv, k_dist = 10, seed = 1) {
  pca_mat <- t(meth_prop)
  rownames(pca_mat) <- cv$sample

  sd_cols <- apply(pca_mat, 2, sd)
  keep <- is.finite(sd_cols) & sd_cols > 0
  pca_mat <- pca_mat[, keep, drop = FALSE]

  set.seed(seed)   # kept for good practice; this fit was previously confirmed fully
                    # deterministic for this p >> n data via a multi-seed check
  pca_fit <- PcaHubert(pca_mat, k = k_dist)
  eig <- pca_fit@eig0

  scores <- as.data.frame(pca_fit@scores)
  colnames(scores) <- paste0("PC", seq_len(ncol(scores)))
  scores$sample <- rownames(scores)

  dat <- merge(as.data.frame(cv), scores, by = "sample", sort = FALSE)
  dat$status <- factor(dat$status)
  dat$colony <- factor(dat$colony)
  treatment_labels <- c("5aza" = "Decitabine", "6aza" = "6-azacytidine", "Control" = "Control")
  dat$treatment <- factor(dat$treatment, levels = names(treatment_labels), labels = treatment_labels)
  dat$status01 <- as.integer(dat$status == "repro")
  dat$stratum  <- interaction(dat$colony, dat$treatment, drop = TRUE)

  list(pca_mat = pca_mat, pca_fit = pca_fit, eig = eig, dat = dat, strata_levels = levels(dat$stratum))
}

## ---- 2.4: standardized centroid distance (Euclidean on standardized PCs = Mahalanobis) ----
add_centroid_distance <- function(d, eig, k) {
  pc_cols <- paste0("PC", seq_len(k))
  sdv <- sqrt(eig[seq_len(k)])
  std <- sweep(as.matrix(d[, pc_cols, drop = FALSE]), 2, sdv, "/")
  ref <- matrix(NA_real_, nrow = nrow(d), ncol = k)
  for (i in seq_len(nrow(d))) {
    same_colony_others <- d$colony == d$colony[i] & d$sample != d$sample[i]
    ref[i, ] <- colMeans(std[same_colony_others, , drop = FALSE])
  }
  colnames(ref) <- paste0("ref_", pc_cols, "_std")
  d$dist_k_ref <- sqrt(rowSums((std - ref)^2))
  d
}

## ---- 2.5: PAIRED relative distance within each stratum (n=2/stratum: one
## sterile, one reproductive). Each individual's distance expressed relative
## to their matched partner - this is what Figure 4/5 panel A actually plots,
## and what the k10_2 draft was missing (see header note, point 5). ----
add_relative_distance <- function(d) {
  d %>%
    group_by(stratum) %>%
    mutate(rel_dist_k_ref = dist_k_ref - rev(dist_k_ref)) %>%
    ungroup() %>%
    as.data.frame()
}

## ---- 2.6: full robustness battery for one label variant's `dat` ----
run_full_battery <- function(dat, strata_levels, label, n_perm = 2000, n_boot = 1000, seed = 1) {

  m_null <- clogit(status01 ~ strata(stratum), data = dat)
  m_dist <- clogit(status01 ~ dist_k_ref + strata(stratum), data = dat)
  a      <- anova(m_null, m_dist)

  lrt_chisq <- a$Chisq[2]
  lrt_p     <- a[[ncol(a)]][2]
  conc      <- summary(m_dist)$concordance[1]

  ## PC1+PC2 linear-model check (see header note, point 6), reusing m_null -
  ## same strata-only null works for both comparisons.
  m_pc12   <- clogit(status01 ~ PC1 + PC2 + strata(stratum), data = dat)
  a_pc12   <- anova(m_null, m_pc12)

  ## Stratified permutation test (2000 reps). Each stratum has exactly one
  ## reproductive and one sterile individual, so permuting within strata =
  ## independently swapping each pair's status labels with probability 0.5.
  set.seed(seed)
  perm_lrt <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    perm_dat <- dat
    for (s in strata_levels) {
      idx <- which(perm_dat$stratum == s)
      if (length(idx) == 2 && sample(c(TRUE, FALSE), 1)) perm_dat$status01[idx] <- rev(perm_dat$status01[idx])
    }
    m0p <- try(clogit(status01 ~ strata(stratum), data = perm_dat), silent = TRUE)
    m1p <- try(clogit(status01 ~ dist_k_ref + strata(stratum), data = perm_dat), silent = TRUE)
    perm_lrt[i] <- if (inherits(m0p, "try-error") || inherits(m1p, "try-error")) NA else anova(m0p, m1p)$Chisq[2]
  }
  perm_n_valid <- sum(!is.na(perm_lrt))
  perm_n_ge    <- sum(perm_lrt >= lrt_chisq, na.rm = TRUE)
  perm_p       <- (perm_n_ge + 1) / (perm_n_valid + 1)

  ## Leave-one-stratum-out
  loso <- lapply(strata_levels, function(s) {
    sub <- dat[dat$stratum != s, ]
    sub$stratum <- droplevels(sub$stratum)
    m0 <- clogit(status01 ~ strata(stratum), data = sub)
    m1 <- clogit(status01 ~ dist_k_ref + strata(stratum), data = sub)
    aa <- anova(m0, m1)
    data.frame(dropped_stratum = s, chisq = aa$Chisq[2], p = aa[[ncol(aa)]][2])
  })
  loso_df <- do.call(rbind, loso)

  ## Stratum-level bootstrap effect-size CI
  set.seed(seed)
  boot_coefs <- sapply(seq_len(n_boot), function(b_i) {
    samp_strata <- sample(strata_levels, length(strata_levels), replace = TRUE)
    boot_dat <- do.call(rbind, lapply(seq_along(samp_strata), function(i) {
      tmp <- dat[dat$stratum == samp_strata[i], ]
      tmp$stratum_boot <- paste0("boot_", i)
      tmp
    }))
    fit <- try(clogit(status01 ~ dist_k_ref + strata(stratum_boot), data = boot_dat), silent = TRUE)
    if (inherits(fit, "try-error")) return(NA_real_)
    coef(fit)["dist_k_ref"]
  })
  boot_ci <- quantile(boot_coefs, c(0.025, 0.975), na.rm = TRUE)

  list(
    label         = label,
    n             = nrow(dat),
    n_strata      = length(strata_levels),
    m_dist        = m_dist,
    lrt_chisq     = lrt_chisq,
    lrt_p         = lrt_p,
    concordance   = conc,
    perm_p        = perm_p,
    perm_n_ge     = perm_n_ge,
    perm_n_valid  = perm_n_valid,
    loso_df       = loso_df,
    loso_range    = range(loso_df$chisq),
    boot_ci       = boot_ci,
    boot_prop_pos = mean(boot_coefs > 0, na.rm = TRUE),
    pc12_chisq    = a_pc12$Chisq[2],
    pc12_p        = a_pc12[[ncol(a_pc12)]][2]
  )
}

# ============================================================
# PART 3: RUN THE k=10 MODEL + FULL BATTERY UNDER ALL THREE LABEL VARIANTS
# Each variant is refit from raw data end to end (its own PCA, its own
# battery) - never sliced or reused across variants.
# ============================================================

variant_results <- list()
variant_fits    <- list()   # keeps matrices/PCA/dat per variant for re-use in
                             # PART 5 (PART 3a), PART 6 (UMAP), PART 8 (sensitivity),
                             # and PART 9 (figures) - all of which run on
                             # PRIMARY_VARIANT only, to avoid recomputing it there.

for (variant in LABEL_VARIANTS) {
  cat(sprintf("\n=== Running k=%d model, label variant = '%s' ===\n", k_dist, variant))

  cv_v   <- build_covar_variant(covar_raw, variant)
  mats_v <- build_meth_matrices(C_raw, N_raw, cv_v)
  fit_v  <- fit_k10_and_build_dat(mats_v$meth_prop, mats_v$cv, k_dist = k_dist)

  fit_v$dat <- add_centroid_distance(fit_v$dat, fit_v$eig, k = k_dist)
  fit_v$dat <- add_relative_distance(fit_v$dat)

  battery_v <- run_full_battery(fit_v$dat, fit_v$strata_levels, label = variant)

  variant_results[[variant]] <- battery_v
  variant_fits[[variant]]    <- c(fit_v, mats_v)

  cat(sprintf(
    "  n=%d, %d strata: chisq=%.3f, p=%.5f, permutation p=%.4f, LOSO %.2f-%.2f, concordance=%.3f\n",
    battery_v$n, battery_v$n_strata, battery_v$lrt_chisq, battery_v$lrt_p, battery_v$perm_p,
    battery_v$loso_range[1], battery_v$loso_range[2], battery_v$concordance
  ))
}

# ============================================================
# PART 4: THREE-WAY LABEL-VARIANT COMPARISON TABLE
# (manuscript Supplementary Table 15)
# ============================================================

label_variant_table <- do.call(rbind, lapply(LABEL_VARIANTS, function(v) {
  r <- variant_results[[v]]
  data.frame(
    variant       = v,
    n             = r$n,
    n_strata      = r$n_strata,
    lrt_chisq     = round(r$lrt_chisq, 2),
    lrt_p         = signif(r$lrt_p, 2),
    permutation_p = signif(r$perm_p, 4),
    loso_min      = round(r$loso_range[1], 1),
    loso_max      = round(r$loso_range[2], 1),
    boot_ci_lower = round(r$boot_ci[1], 2),
    boot_ci_upper = round(r$boot_ci[2], 2),
    concordance   = round(r$concordance, 3)
  )
}))
cat("\n============ THREE-WAY LABEL-VARIANT COMPARISON ============\n")
print(label_variant_table)
write.csv(label_variant_table, file.path(out_dir, "label_variant_comparison_omnibus.csv"), row.names = FALSE)

# ============================================================
# PART 5: [OPTIONAL] EIGENVALUE SPECTRUM + PER-COMPONENT (PC1-PC10)
# PERMUTATION NULLS, ROBUST PCA
# Ported back from the original univariate script's PART 3a (dropped
# entirely in _centroid_k10_2.R - see header note, point 3), then extended
# from a PC1-only test to a per-component test across all k_dist components
# (see header note, point 7): a proper parallel-analysis-style procedure,
# where each permutation replicate is fit once at the full k=10 and every
# component's observed variance share is compared against its own null
# distribution, rather than testing PC1 alone. Runs on the primary
# variant's already-fitted k=10 PCA; does NOT refit a separate, larger PCA
# the way the original script did (that script requested a bigger fit here
# than the model itself used - the "k_extract must equal the number of
# components used, exactly" lesson applies to this block too).
# ============================================================

if (RUN_PART3A) {

  primary_fit     <- variant_fits[[PRIMARY_VARIANT]]
  pca_mat_primary <- primary_fit$pca_mat
  eig_full        <- primary_fit$eig   # eig0 from the k=10 fit itself

  pc_variance_table <- data.frame(
    PC                  = paste0("PC", seq_along(eig_full)),
    eigenvalue          = eig_full,
    proportion_variance = eig_full / sum(eig_full),
    cumulative_variance = cumsum(eig_full / sum(eig_full))
  )
  cat("\n============ PC VARIANCE EXPLAINED (k=10 fit, primary variant) ============\n")
  print(pc_variance_table)
  write.csv(pc_variance_table, file.path(out_dir, "pc_variance_explained_k10_omnibus.csv"), row.names = FALSE)

  ## ---- Robust PCA permutation null, ALL k_dist COMPONENTS (not just PC1) ----
  ## Still only n_perm_robust total PcaHubert calls, same as the old PC1-only
  ## version - each replicate is fit once at k=k_dist and its full eigenvalue
  ## vector is kept, rather than fitting once per component.
  ##
  ## NOTE: @eig0 (both here and for the primary fit's eig_full above) returns
  ## the FULL diagnostic eigenvalue spectrum from PcaHubert's internal
  ## dimension-reduction step, NOT just the k_dist components actually
  ## retained/requested - this is why pc_variance_table above can run past
  ## PC10. Every use of eig0 below is written to cope with that (indexing/
  ## subsetting to the first k_dist entries where a fixed-length result is
  ## needed), rather than assuming length(eig0) == k_dist.
  set.seed(1)
  n_perm_robust    <- 49
  prop_perm_robust <- matrix(NA_real_, nrow = n_perm_robust, ncol = k_dist)
  for (i in seq_len(n_perm_robust)) {
    perm_mat <- apply(pca_mat_primary, 2, sample)
    pc_perm  <- PcaHubert(perm_mat, k = k_dist)
    prop_perm_robust[i, ] <- (pc_perm@eig0 / sum(pc_perm@eig0))[seq_len(k_dist)]
  }

  obs_prop_robust <- eig_full / sum(eig_full)
  p_perm_by_pc <- sapply(seq_len(k_dist), function(j) {
    (sum(prop_perm_robust[, j] >= obs_prop_robust[j]) + 1) / (n_perm_robust + 1)
  })

  pc_permutation_table <- data.frame(
    PC                  = paste0("PC", seq_len(k_dist)),
    proportion_variance = round(obs_prop_robust[seq_len(k_dist)], 4),   # eig_full/obs_prop_robust
    permutation_p       = signif(p_perm_by_pc, 3)                      # can be LONGER than k_dist -
  )                                                                     # @eig0 is the full diagnostic
                                                                         # spectrum, not just the k
                                                                         # retained components (see
                                                                         # pc_variance_table above,
                                                                         # which already showed this)
  cat(sprintf("\n============ ROBUST PCA PER-COMPONENT PERMUTATION NULL (n=%d reps, k=%d) ============\n",
              n_perm_robust, k_dist))
  print(pc_permutation_table)
  write.csv(pc_permutation_table, file.path(out_dir, "pc_permutation_null_k10_omnibus.csv"), row.names = FALSE)

  cat(sprintf("\nComponents with permutation p < 0.05: %s\n",
              paste(pc_permutation_table$PC[p_perm_by_pc < 0.05], collapse = ", ")))
  cat(sprintf("Cumulative variance, PC1-%d: %.1f%%\n", k_dist, 100 * pc_variance_table$cumulative_variance[k_dist]))

  ## ---- [OPTIONAL] Classical (non-robust) PCA permutation null, PC1 only ----
  ## Off by default (RUN_CLASSICAL_PERM_NULL) - see header note, point 7, for
  ## why. Code kept as-is, just gated, in case you want it back.
  if (RUN_CLASSICAL_PERM_NULL) {
    set.seed(1)
    n_perm_classical <- 199
    pc1_prop_perm_classical <- numeric(n_perm_classical)
    for (i in seq_len(n_perm_classical)) {
      perm_mat <- apply(pca_mat_primary, 2, sample)
      ev_perm  <- prcomp(perm_mat, center = TRUE, scale. = FALSE)$sdev^2
      pc1_prop_perm_classical[i] <- ev_perm[1] / sum(ev_perm)
    }
    ev_obs_classical       <- prcomp(pca_mat_primary, center = TRUE, scale. = FALSE)$sdev^2
    pc1_prop_obs_classical <- ev_obs_classical[1] / sum(ev_obs_classical)
    p_perm_classical <- (sum(pc1_prop_perm_classical >= pc1_prop_obs_classical) + 1) / (n_perm_classical + 1)

    cat(sprintf(
      "\n[Classical PCA cross-check] PC1 variance share = %.1f%% (perm p %s)\n",
      100 * pc1_prop_obs_classical,
      ifelse(p_perm_classical < 0.005, "< 0.005", sprintf("= %.3f", p_perm_classical))
    ))
  } else {
    cat("\n[Classical PCA permutation cross-check skipped: RUN_CLASSICAL_PERM_NULL = FALSE]\n")
  }

} else {
  cat("\n[PART 5 skipped: RUN_PART3A = FALSE]\n")
}

# ============================================================
# PART 6: UMAP (k=10) - colony-clustering figure, primary variant
# Ported back from the original script's PART 4b (dropped entirely in
# _centroid_k10_2.R - see header note, point 4). Now naturally uses all 10
# PCs from the k=10 fit; the original script's UMAP was capped at 5 PCs as
# a side-effect of that script's k_extract=5 PCA, not a deliberate choice.
# ============================================================

primary_dat <- variant_fits[[PRIMARY_VARIANT]]$dat
primary_pca <- variant_fits[[PRIMARY_VARIANT]]$pca_fit

n_umap_pcs <- min(k_dist, ncol(primary_pca@scores))   # = 10
pc_umap <- scale(primary_pca@scores[, 1:n_umap_pcs, drop = FALSE])

set.seed(1)
um <- umap(pc_umap, n_neighbors = 6, min_dist = 0.3, metric = "euclidean")

primary_dat$UMAP1 <- um[, 1]
primary_dat$UMAP2 <- um[, 2]
primary_dat$fill_col <- ifelse(primary_dat$status == "repro", as.character(primary_dat$colony), NA)
primary_dat$fill_col <- factor(primary_dat$fill_col)

p_global <- ggplot(primary_dat, aes(UMAP1, UMAP2, group = interaction(colony, treatment))) +
  geom_line(aes(linetype = treatment), colour = "grey30", linewidth = 0.7, alpha = 0.7) +
  geom_point(aes(colour = colony, fill = fill_col), shape = 21, size = 1.8, alpha = 0.9, stroke = 0.8) +
  scale_fill_discrete(na.value = "transparent", guide = "none") +
  theme_classic(base_size = 11) +
  labs(x = "UMAP1", y = "UMAP2", colour = "Colony", linetype = "Treatment")

p_facet <- ggplot(primary_dat, aes(UMAP1, UMAP2, group = interaction(colony, treatment))) +
  geom_line(aes(linetype = treatment), colour = "grey30", linewidth = 0.7, alpha = 0.7) +
  geom_point(aes(colour = colony, fill = fill_col), shape = 21, size = 1.8, alpha = 0.9, stroke = 0.8) +
  facet_wrap(~ colony, scales = "free") +
  scale_colour_discrete(guide = "none") +
  scale_fill_discrete(na.value = "transparent", guide = "none") +
  theme_classic(base_size = 11) +
  theme(strip.text = element_blank()) +
  labs(x = "UMAP1", y = "UMAP2", linetype = "Treatment") +
  guides(linetype = "none") +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank()
  )

combined_umap_plot <- (p_global + p_facet) +
  plot_layout(widths = c(1.2, 1), guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.tag = element_text(size = 14, face = "bold")
  )
ggsave(file.path(out_dir, "facet_umap_k10_omnibus.pdf"), combined_umap_plot, width = 11.4, height = 4.9)

## Supplementary common-axis version (fixed scales across colony panels)
p_facet_common <- ggplot(primary_dat, aes(UMAP1, UMAP2, group = interaction(colony, treatment))) +
  geom_line(aes(linetype = treatment), colour = "grey30", linewidth = 0.7, alpha = 0.7) +
  geom_point(aes(colour = colony, fill = fill_col), shape = 21, size = 1.8, alpha = 0.9, stroke = 0.8) +
  facet_wrap(~ colony) +   # default scales = "fixed"
  scale_colour_discrete(guide = "none") +
  scale_fill_discrete(na.value = "transparent", guide = "none") +
  theme_classic(base_size = 11) +
  labs(x = "UMAP1", y = "UMAP2", linetype = "Treatment") +
  guides(linetype = "none")
ggsave(file.path(out_dir, "facet_umap_common_axes_k10_omnibus_supp.pdf"), p_facet_common, width = 9, height = 6.5)

# ============================================================
# PART 7: PC1+PC2 LINEAR-MODEL CHECK, ACROSS ALL THREE VARIANTS
# Recomputed fresh from each variant's own k=10 PCA fit (see header note,
# point 6) - do not carry forward the old manuscript's 0.74/0.70, which was
# computed under a k_extract=5 fit.
# ============================================================

pc12_table <- do.call(rbind, lapply(LABEL_VARIANTS, function(v) {
  r <- variant_results[[v]]
  data.frame(variant = v, chisq = round(r$pc12_chisq, 2), p = round(r$pc12_p, 2))
}))
cat("\n============ PC1+PC2 LINEAR MODEL vs. STRATA-ONLY NULL ============\n")
print(pc12_table)
write.csv(pc12_table, file.path(out_dir, "pc1pc2_linear_check_omnibus.csv"), row.names = FALSE)

# ============================================================
# PART 8: [OPTIONAL] 5-VARIANT SENSITIVITY SWEEP
# (transform / PCA method / imputation), primary variant only
# (manuscript Supplementary Table 16).
# ============================================================

run_centroid_pipeline <- function(prop_mat_samples_x_cpg, cv, use_robust = TRUE,
                                   k_dist_local = 10, label = "") {
  sd_cols <- apply(prop_mat_samples_x_cpg, 2, sd)
  keep <- is.finite(sd_cols) & sd_cols > 0
  X <- prop_mat_samples_x_cpg[, keep, drop = FALSE]
  k <- min(k_dist_local, nrow(X) - 1, ncol(X))   # request EXACTLY k_dist_local components

  if (use_robust) {
    pc_fit <- PcaHubert(X, k = k)
    scores <- as.data.frame(pc_fit@scores)
    eig <- pc_fit@eig0
  } else {
    pc_fit <- prcomp(X, center = TRUE, scale. = FALSE)
    scores <- as.data.frame(pc_fit$x[, 1:k, drop = FALSE])
    eig <- pc_fit$sdev^2
  }
  colnames(scores) <- paste0("PC", seq_len(ncol(scores)))
  scores$sample <- rownames(X)

  d <- merge(as.data.frame(cv), scores, by = "sample", sort = FALSE)
  d$status01  <- as.integer(factor(d$status) == "repro")
  d$colony    <- factor(d$colony)
  d$treatment <- factor(d$treatment)
  d$stratum   <- interaction(d$colony, d$treatment, drop = TRUE)
  d <- add_centroid_distance(d, eig, k = k_dist_local)

  m0 <- clogit(status01 ~ strata(stratum), data = d)
  m1 <- clogit(status01 ~ dist_k_ref + strata(stratum), data = d)
  a  <- anova(m0, m1)

  list(label = label, n_cpg = ncol(X), chisq = a$Chisq[2], p = a[[ncol(a)]][2])
}

if (RUN_SENSITIVITY) {

  pm <- variant_fits[[PRIMARY_VARIANT]]   # primary variant's cv + matrices

  asin_mat <- asin(sqrt(pmin(pmax(pm$meth_prop, 0), 1)))

  elogit_mat <- log((pm$C_mat + 0.5) / (pm$N_mat - pm$C_mat + 0.5))
  elogit_mat[pm$N_mat == 0] <- NA
  if (anyNA(elogit_mat)) {
    rm_e <- rowMeans(elogit_mat, na.rm = TRUE)
    na_rows_e <- which(rowSums(is.na(elogit_mat)) > 0)
    for (i in na_rows_e) elogit_mat[i, is.na(elogit_mat[i, ])] <- rm_e[i]
  }

  keep_cc <- rowSums(is.na(pm$meth_prop_raw)) == 0
  complete_case_mat <- pm$meth_prop_raw[keep_cc, , drop = FALSE]
  cat(sprintf("Complete-case filter: %d of %d CpGs retained (primary variant)\n",
              nrow(complete_case_mat), nrow(pm$meth_prop_raw)))

  results_centroid_sens <- list(
    original_robust      = run_centroid_pipeline(t(pm$meth_prop), pm$cv, TRUE,
                              label = "Original: proportions, mean-imputed, robust PCA"),
    asin_robust          = run_centroid_pipeline(t(asin_mat), pm$cv, TRUE,
                              label = "Arcsine-sqrt transform, robust PCA"),
    elogit_robust        = run_centroid_pipeline(t(elogit_mat), pm$cv, TRUE,
                              label = "Empirical logit transform, robust PCA"),
    original_classical   = run_centroid_pipeline(t(pm$meth_prop), pm$cv, FALSE,
                              label = "Original proportions, CLASSICAL (non-robust) PCA"),
    complete_case_robust = run_centroid_pipeline(t(complete_case_mat), pm$cv, TRUE,
                              label = "Complete-case (no imputation), robust PCA")
  )

  centroid_sens_summary <- do.call(rbind, lapply(results_centroid_sens, function(r) {
    data.frame(label = r$label, n_cpg = r$n_cpg, chisq = round(r$chisq, 3), p = round(r$p, 4))
  }))
  rownames(centroid_sens_summary) <- NULL
  cat("\n============ SENSITIVITY SWEEP (k=10 centroid distance, primary variant) ============\n")
  print(centroid_sens_summary)
  write.csv(centroid_sens_summary,
            file.path(out_dir, "distance_model_sensitivity_summary_centroid_k10_omnibus.csv"),
            row.names = FALSE)

} else {
  cat("\n[PART 8 skipped: RUN_SENSITIVITY = FALSE]\n")
}

# ============================================================
# PART 9: FIGURES (primary variant)
# Predicted-probability curve now uses the PAIRED relative distance
# (rel_dist_k_ref), matching Figure 4/5 panel A as described in the
# manuscript - see header note, point 5.
# ============================================================

primary_dat2 <- variant_fits[[PRIMARY_VARIANT]]$dat   # already has dist_k_ref + rel_dist_k_ref
b_dist <- coef(variant_results[[PRIMARY_VARIANT]]$m_dist)

set.seed(42)
primary_dat2$jitter_status <- jitter(primary_dat2$status01, amount = 0.02)

rel_seq <- seq(min(primary_dat2$rel_dist_k_ref), max(primary_dat2$rel_dist_k_ref), length.out = 200)
pred_grid_rel <- data.frame(rel_dist_k_ref = rel_seq)
pred_grid_rel$prob <- plogis(b_dist["dist_k_ref"] * pred_grid_rel$rel_dist_k_ref)

p_dist_curve_rel <- ggplot() +
  geom_point(data = primary_dat2, aes(x = rel_dist_k_ref, y = jitter_status, colour = colony),
             size = 2.2, alpha = 0.8) +
  geom_line(data = pred_grid_rel, aes(x = rel_dist_k_ref, y = prob), colour = "black", linewidth = 1.1) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  theme_classic(base_size = 12) +
  labs(x = "Distance from colony reference, relative to pair partner",
       y = "Reproductive status / predicted probability", colour = "Colony") +
  guides(colour = "none")

set.seed(42)   # reproducible jitter, matching the seed already used above for jitter_status
primary_dat2$x_jit <- as.numeric(primary_dat2$status) +
  runif(nrow(primary_dat2), -0.08, 0.08)

p_dist_box <- ggplot(primary_dat2, aes(x = status, y = dist_k_ref, fill = status)) +
  geom_boxplot(width = 0.5, alpha = 0.6, outlier.shape = NA) +
  geom_line(aes(x = x_jit, group = stratum), colour = "grey60", linewidth = 0.4, alpha = 0.6) +
  geom_point(aes(x = x_jit, colour = colony), size = 2.2, alpha = 0.8) +
  scale_fill_manual(values = c("repro" = "black", "sterile" = "white")) +
  theme_classic(base_size = 12) +
  labs(x = "Reproductive status", y = "k=10 centroid distance", colour = "Colony") +
  guides(fill = "none")

combined_dist_plot <- (p_dist_curve_rel + p_dist_box) +
  plot_layout(widths = c(1.15, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 14, face = "bold"))
ggsave(file.path(out_dir, "pc1_dist_combined_centroid_k10_omnibus.pdf"), combined_dist_plot, width = 12, height = 5.5)
ggsave(file.path(out_dir, "pc1_dist_boxplot_centroid_k10_omnibus.pdf"), p_dist_box, width = 6, height = 5)

# ============================================================
# PART 10: OMNIBUS VERDICT
# ============================================================

pv <- variant_results[[PRIMARY_VARIANT]]
cat(
  "\n================ k=10 CENTROID-DISTANCE MODEL: OMNIBUS VERDICT ================\n",
  sprintf("Primary variant ('%s'): n = %d, %d strata\n", PRIMARY_VARIANT, pv$n, pv$n_strata),
  sprintf("LRT (vs. no-covariate null): chisq = %.3f, p = %.5f\n", pv$lrt_chisq, pv$lrt_p),
  sprintf("Permutation LRT:             p = %.4f (%d/%d permutations met or exceeded observed)\n",
          pv$perm_p, pv$perm_n_ge, pv$perm_n_valid),
  sprintf("LOSO chisq range (%d):        %.2f to %.2f\n", pv$n_strata, pv$loso_range[1], pv$loso_range[2]),
  sprintf("Bootstrap effect-size 95%% CI: %.3f to %.3f (%.0f%% of replicates positive)\n",
          pv$boot_ci[1], pv$boot_ci[2], 100 * pv$boot_prop_pos),
  sprintf("Concordance:                  %.3f\n", pv$concordance),
  sprintf("PC1+PC2 linear check:         chisq = %.2f, p = %.2f (recomputed at k=10 - see header note 6)\n",
          pv$pc12_chisq, pv$pc12_p),
  "\nThree-way label-variant comparison: see label_variant_table above / label_variant_comparison_omnibus.csv\n",
  if (RUN_PART3A) "Eigenvalue spectrum / per-component permutation nulls: see PART 5 output above.\n" else "PART 5 (eigenvalue spectrum) was skipped this run.\n",
  if (RUN_SENSITIVITY) "5-variant sensitivity sweep: see centroid_sens_summary above.\n" else "PART 8 (sensitivity sweep) was skipped this run.\n",
  "===================================================================================\n\n"
)


# ------------------- GLOBAL TREATMENT EFFECT ON METHYLATION STRUCTURE ----------
# PERMANOVA on robust PC space, blocked by colony
library(vegan)

pf <- variant_fits[[PRIMARY_VARIANT]]
K  <- min(k_dist, ncol(pf$pca_fit@scores))
D  <- dist(scale(pf$pca_fit@scores[, 1:K, drop = FALSE]))
adonis2(D ~ treatment, data = pf$dat, permutations = 999, strata = pf$dat$colony)

# ============================================================
# END
# ============================================================
