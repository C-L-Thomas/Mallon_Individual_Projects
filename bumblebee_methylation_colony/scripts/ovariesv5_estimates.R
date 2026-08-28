## =========================================================
## Ovary analysis: status, score, and length
## v5 -- adds pairwise treatment contrasts (absolute estimates
## + 95% CIs) for Reviewer 1, major comment 8:
## "Present treatment effects with estimates, uncertainty and
## appropriate caution about non-significant tests."
##
## Everything above the "PAIRWISE CONTRASTS" blocks is unchanged
## from ovariesv4.R. New code is marked NEW.
## =========================================================

library(lme4)
library(ordinal)
library(censReg)
library(car)
library(emmeans)
library(dplyr)
library(ggplot2)
library(patchwork)

## ---- Data prep ----

# Work vs home machines may use different OneDrive account folders.
find_existing <- function(paths) {
  hit <- paths[dir.exists(paths)]
  if (length(hit) == 0) {
    stop("None of the candidate base directories exist:\n", paste(paths, collapse = "\n"))
  }
  hit[1]
}
base_dir <- find_existing(c(
  "/Users/emb3/Library/CloudStorage/OneDrive-UniversityofLeicester/projects/Eamonn_Bumblebee/copy_for_paper/revision",
  "/Users/ebm3/Library/CloudStorage/OneDrive-UniversityofLeicester/projects/Eamonn_Bumblebee/copy_for_paper/revision"
))
data_path   <- file.path(base_dir, "data", "Reproductive_Data_christian.csv")
output_dir  <- file.path(base_dir, "outputs")

dat <- read.csv(
  data_path,
  skip = 1,               # first line is a title row ("Supplemental Information Table 2"), not data
  stringsAsFactors = TRUE
)

lvl <- c("Control", "6-azacytidine", "Decitabine")
dat$Drug   <- factor(dat$Drug, levels = lvl)
dat$Colony <- factor(dat$Colony)
dat$Box    <- factor(dat$Box)

dat$y <- as.integer(dat$Reproductive == "Reproductive")
dat$Score_ord <- factor(dat$Score, levels = 0:4, ordered = TRUE)

treat_cols <- c(
  "Control"       = "#333333",
  "6-azacytidine" = "#1b9e77",
  "Decitabine"    = "#d95f02"
)

## =========================================================
## 1. Reproductive status (binary) -- GLMM
## =========================================================

m_status <- glmer(
  y ~ Drug + (1 | Colony/Box),
  family  = binomial,
  data    = dat,
  control = glmerControl(optimizer = "bobyqa")
)

summary(m_status)
car::Anova(m_status, type = "II")

emm_status <- emmeans(m_status, ~ Drug, type = "response")
pairs(emm_status, adjust = "tukey")

status_df <- as.data.frame(emm_status) %>%
  transmute(Drug, estimate = prob, lower = asymp.LCL, upper = asymp.UCL) %>%
  mutate(Drug = factor(Drug, levels = lvl))

## ---- NEW: pairwise contrasts, two scales ----
## (a) Odds ratio scale -- what pairs() gives by default on a
##     type="response" grid from a logit-link model.
status_pairs_OR <- pairs(emm_status, adjust = "tukey")
confint(status_pairs_OR)

## (b) Absolute risk difference (probability scale) -- what the
##     reviewer means by "absolute model estimates". regrid()
##     drops the response transformation tag so pairs() takes
##     plain differences of probabilities instead of ratios.
status_pairs_RD <- emm_status %>%
  regrid() %>%
  pairs(adjust = "tukey")
confint(status_pairs_RD)

## Column names to read off: `estimate`/`odds.ratio` (whichever
## appears), `asymp.LCL`, `asymp.UCL`, `p.value`. Run
## as.data.frame(status_pairs_RD) if unsure.

## =========================================================
## 2. Reproductive score (0-4 ordinal) -- CLMM
## =========================================================

m_score <- clmm(
  Score_ord ~ Drug + (1 | Colony/Box),
  data = dat,
  link = "logit"
)

summary(m_score)

m_score_null <- clmm(
  Score_ord ~ 1 + (1 | Colony/Box),
  data = dat,
  link = "logit"
)
anova(m_score_null, m_score)

emm_score <- emmeans(m_score, ~ Drug, mode = "latent")
pairs(emm_score, adjust = "tukey")

score_df <- as.data.frame(emm_score) %>%
  transmute(Drug, estimate = emmean, lower = asymp.LCL, upper = asymp.UCL) %>%
  mutate(Drug = factor(Drug, levels = lvl))


emm_cat  <- emmeans(m_score, ~ Drug, mode = "mean.class")
score_df <- as.data.frame(emm_cat) %>%
  transmute(Drug,
            estimate = mean.class,
            lower    = asymp.LCL,
            upper    = asymp.UCL) %>%
  mutate(Drug = factor(Drug, levels = lvl))

## ---- NEW: pairwise contrasts, two scales ----
## (a) Absolute scale: difference in expected score (0-4 units).
##     This is the scale to quote as "absolute model estimates".
score_pairs_absolute <- pairs(emm_cat, adjust = "tukey")
confint(score_pairs_absolute)

## (b) Latent (log-cumulative-odds) scale, for readers who want
##     the proportional-odds-style effect size.
score_pairs_latent <- pairs(emm_score, adjust = "tukey")
confint(score_pairs_latent)

## =========================================================
## 3. Ovary length -- Tobit (left-censored at 0)
## =========================================================

m_length <- censReg(
  Length ~ Drug + Colony,
  left = 0,
  data = dat
)

summary(m_length)
margEff(m_length)

## Extract Drug estimates with 95% CIs on the latent scale.
## censReg doesn't play with emmeans, so we pull coefficients directly
## and compute marginal means relative to the Control reference.
coefs  <- coef(summary(m_length))
vc     <- vcov(m_length)

# Intercept = Control mean (at reference Colony level)
mu_control  <- coefs["(Intercept)", "Estimate"]
se_control  <- coefs["(Intercept)", "Std. error"]

# Other treatments = intercept + Drug coefficient
b_aza  <- coefs["Drug6-azacytidine", "Estimate"]
b_dec  <- coefs["DrugDecitabine", "Estimate"]

mu_aza <- mu_control + b_aza
mu_dec <- mu_control + b_dec

# SEs using the covariance matrix
se_aza <- sqrt(vc["(Intercept)", "(Intercept)"] +
               vc["Drug6-azacytidine", "Drug6-azacytidine"] +
               2 * vc["(Intercept)", "Drug6-azacytidine"])
se_dec <- sqrt(vc["(Intercept)", "(Intercept)"] +
               vc["DrugDecitabine", "DrugDecitabine"] +
               2 * vc["(Intercept)", "DrugDecitabine"])

length_df <- data.frame(
  Drug     = factor(lvl, levels = lvl),
  estimate = c(mu_control, mu_aza, mu_dec),
  se       = c(se_control, se_aza, se_dec)
) %>%
  mutate(
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se
  )

## ---- NEW: pairwise treatment contrasts with 95% CIs ----
## b_aza and b_dec ARE the Control-referenced contrasts already
## (Control is the reference level of Drug), so their SEs come
## straight off the coefficient table -- no need to combine with
## the intercept variance the way length_df does for the absolute
## group means above. The third contrast (Decitabine vs
## 6-azacytidine) has to be built by hand from the covariance
## matrix since censReg has no emmeans method.

se_aza_vs_control <- sqrt(vc["Drug6-azacytidine", "Drug6-azacytidine"])
se_dec_vs_control <- sqrt(vc["DrugDecitabine", "DrugDecitabine"])

cov_aza_dec <- vc["Drug6-azacytidine", "DrugDecitabine"]
se_dec_vs_aza <- sqrt(vc["DrugDecitabine", "DrugDecitabine"] +
                      vc["Drug6-azacytidine", "Drug6-azacytidine"] -
                      2 * cov_aza_dec)

## IMPORTANT: use the same Tukey correction here as for the status
## and score models above (status_pairs_RD, score_pairs_absolute both
## used adjust = "tukey"). A plain z-based CI on these three contrasts
## is NOT comparable to the other two models and will make one
## contrast look "significant" when it wouldn't survive the same
## correction applied elsewhere in the paper.
tukey_mult <- qtukey(0.95, nmeans = 3, df = Inf) / sqrt(2)  # ~2.344, matches the
                                                             # emmeans Tukey CIs above

length_pairs <- data.frame(
  contrast = c("6-azacytidine - Control",
               "Decitabine - Control",
               "Decitabine - 6-azacytidine"),
  estimate = c(b_aza, b_dec, b_dec - b_aza),
  se       = c(se_aza_vs_control, se_dec_vs_control, se_dec_vs_aza)
) %>%
  mutate(
    # unadjusted, for reference only -- do not report these alongside
    # Tukey-adjusted CIs from the other two models
    lower_unadj = estimate - 1.96 * se,
    upper_unadj = estimate + 1.96 * se,
    z           = estimate / se,
    p_unadj     = 2 * pnorm(-abs(z)),
    # Tukey-adjusted -- use these in the manuscript
    lower_tukey = estimate - tukey_mult * se,
    upper_tukey = estimate + tukey_mult * se
  )

length_pairs
## Estimates are on the latent (uncensored) Gaussian scale, i.e.
## mm of terminal oocyte length as if measurement were not
## censored at 0. Same scale as mu_control/mu_aza/mu_dec above,
## so it's directly comparable to Figure 2C.

library(lmtest)

m_length      <- censReg(Length ~ Drug + Colony, left = 0, data = dat)
m_length_null <- censReg(Length ~ Colony,        left = 0, data = dat)

lrtest(m_length_null, m_length)

## =========================================================
## Console summary of everything the reviewer asked for
## =========================================================
cat("\n--- Reproductive status: absolute risk differences ---\n")
print(confint(status_pairs_RD))

cat("\n--- Reproductive score: differences in expected score (0-4) ---\n")
print(confint(score_pairs_absolute))

cat("\n--- Terminal oocyte length: differences on latent (mm) scale ---\n")
print(length_pairs)

## =========================================================
## Unified figure: three panels, one measure of ovary per panel
## =========================================================

base_theme <- theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.title = element_text(size = 10, face = "plain")
  )

p_status <- ggplot(status_df, aes(Drug, estimate, colour = Drug)) +
  geom_point(size = 2.8) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15) +
  scale_colour_manual(values = treat_cols) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = NULL, y = "P(reproductive)",
       title = "Reproductive status (GLMM)") +
  base_theme



p_score <- ggplot(score_df, aes(Drug, estimate, colour = Drug)) +
  geom_jitter(
    data = dat,
    aes(x = Drug, y = as.numeric(as.character(Score)), colour = Drug),
    width = 0.08, height = 0.08, alpha = 0.25, inherit.aes = FALSE
  ) +
  geom_point(size = 2.8) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15) +
  scale_colour_manual(values = treat_cols) +
  scale_y_continuous(
    breaks = 0:4,
    limits = c(-0.3, 4.3)
  ) +
  labs(x = NULL, y = "Reproductive score (0–4)",
       title = "Reproductive score (CLMM)") +
  base_theme





p_length <- ggplot(length_df, aes(Drug, estimate, colour = Drug)) +
  geom_jitter(
    data = dat,
    aes(x = Drug, y = Length, colour = Drug),
    width = 0.08, height = 0, alpha = 0.25, inherit.aes = FALSE
  ) +
  geom_point(size = 2.8) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15) +
  scale_colour_manual(values = treat_cols) +
  labs(x = NULL, y = "Ovary length (mm)",
       title = "Ovary length (Tobit)") +
  base_theme

combined <- p_status + p_score + p_length +
  plot_annotation(tag_levels = "A")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(
  file.path(output_dir, "ovary_three_panel.pdf"),
  combined,
  width  = 10,
  height = 4
)
