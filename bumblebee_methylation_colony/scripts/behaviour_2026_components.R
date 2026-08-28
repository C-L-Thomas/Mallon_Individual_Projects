
library(ggplot2)
library(glmmTMB)
library(dplyr)
library(splines)
library(emmeans)




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

data_dir   <- file.path(base_dir, "data")
output_dir <- file.path(base_dir, "outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## Aggression_Data.csv (unlike complied_aggression.csv) has a title row
## ("Supplemental Information Table 1") above the real header, hence skip = 1.
comp_data <- read.csv(file.path(data_dir, "Aggression_Data.csv"), skip = 1)

comp_data$Treatment   <- factor(comp_data$Treatment)
comp_data$Colony <- factor(comp_data$Colony)
comp_data$Box.Name    <- factor(comp_data$Box.Name)
comp_data$Batch <- factor(comp_data$Batch)
comp_data$Date <- as.Date(comp_data$Date, format = "%d/%m/%Y")

## NOTE: Aggression_Data.csv's Date column is not consistently formatted
## (some rows DD/MM/YYYY, some MM/DD/YYYY -- confirmed by cross-checking
## against complied_aggression.csv, where the same rows are consistent
## DD/MM/YYYY). The parse above will silently misread the ambiguous rows.
## This does not affect anything below: DayNum (next line) comes from the
## "D1".."D6" text labels in the Day column, not from Date, and Date/Day
## (recomputed below) are only used for the summary() sanity check, not
## the models or plots. If Date is ever needed for a real calendar-date
## calculation, fix the source column first rather than trusting this parse.

comp_data$DayNum <- as.integer(sub("D", "", comp_data$Day))

comp_data$Day <- as.integer(comp_data$Date - min(comp_data$Date, na.rm = TRUE))
summary(comp_data$Day)

comp_data$Treatment <- relevel(comp_data$Treatment, ref = "Control")

## CHECK: confirm the composite Total column and the three component
## columns below actually match your CSV headers before running.
print(names(comp_data))

m_total_spline <- glmmTMB(
  Total ~ Treatment + Time + ns(DayNum, df = 3) +
    (1 | Colony/Box.Name),
  data = comp_data,
  family = nbinom2
)
summary(m_total_spline) #Aggression changes strongly and non-linearly across experimental days.
performance::check_overdispersion(m_total_spline)

#Aggressive behaviour counts were analysed using negative binomial mixed-effects models with treatment, time of day, and a natural spline for experimental day as fixed effects, and box nested within colony as random intercepts to account for repeated measures.

#Aggressive behaviour differed significantly among treatments. Relative to control boxes, aggression increased by approximately 36% under 6-azacytidine and more than doubled under Decitabine (negative binomial mixed-effects model, both p < 0.05). Aggression also varied strongly and non-linearly across experimental days, consistent with habituation effects, but did not differ systematically across times of day.

emm_treat <- emmeans(
  m_total_spline,
  ~ Treatment,
  type = "response"
)

emm_treat

## reverse = TRUE flips ratio/estimate direction so contrasts read as
## treatment/reference (e.g. Decitabine/Control) instead of emmeans'
## default reference/treatment ordering -- matches how results are
## written up in prose throughout this script.
pairs(emm_treat, reverse = TRUE)

#Model-estimated aggressive behaviour increased from a mean of 1.48 events per observation in control boxes to 2.02 under 6-azacytidine and 3.28 under Decitabine. Pairwise comparisons indicated a strong increase under Decitabine relative to both control (rate ratio ~ 2.21, p < 0.001) and 6-azacytidine (rate ratio ~ 1.62, p < 0.001). The increase under 6-azacytidine relative to control was more modest (rate ratio ~ 1.36) and marginal after Tukey correction (p = 0.067).


## =========================================================
## NEW: component behaviours (Reviewer 1)
##
## "The composite total of attack, darting and humming is useful, but
## the component behaviours ... should be checked and reported, at
## least in the supplement, because they may not be biologically
## equivalent."
##
## The composite Total is an unweighted sum of three behaviours that
## differ a lot biologically: attacking is contact aggression, darting
## and humming are non-contact threat displays. A treatment could
## plausibly shift one component without shifting the others, and the
## sum would hide that. This section fits the identical model
## structure used for Total to each component separately, so any
## divergence between components is directly visible rather than
## averaged away.
## =========================================================

## Edit these three names to match the actual column headers in
## Aggression_Data.csv (checked via names(comp_data) above).
component_vars <- c("Attack", "Darting", "Humming")

stopifnot(all(component_vars %in% names(comp_data)))

fit_component_model <- function(response_var, data = comp_data) {
  f <- as.formula(
    paste0(response_var, " ~ Treatment + Time + ns(DayNum, df = 3) + (1 | Colony/Box.Name)")
  )
  glmmTMB(f, data = data, family = nbinom2)
}

component_models <- lapply(component_vars, fit_component_model)
names(component_models) <- component_vars

## Same diagnostics you already run on m_total_spline, once per
## component. Watch for convergence warnings and overdispersion flags
## in particular -- a rarer behaviour (e.g. humming, if it has a lot
## of zero counts) may not tolerate the same nbinom2 spec as Total and
## could need family = "compois" or a simpler day term.
lapply(component_models, summary)
lapply(component_models, performance::check_overdispersion)

## Marginal means (response/count scale) and Tukey pairwise contrasts
## per component -- same scale and adjustment as the Total model
## above, so they're directly comparable.
component_emm   <- lapply(component_models, function(m) emmeans(m, ~ Treatment, type = "response"))
component_pairs <- lapply(component_emm, pairs, adjust = "tukey", reverse = TRUE)

## Same emmeans/pairs block as the Total model, written out per
## component so the console output (with exact Tukey-adjusted
## p-values) matches the format of the existing Total narrative
## comment above and can be turned into matching sentences directly.

emm_attack <- emmeans(component_models[["Attack"]], ~ Treatment, type = "response")
emm_attack
pairs(emm_attack, reverse = TRUE)
#Model-estimated attacking was 0.25 events per observation in control boxes, 0.23 under
#6-azacytidine, and 0.54 under Decitabine. Decitabine vs control: rate ratio ~ 2.17,
#p = [fill in]. Decitabine vs 6-azacytidine: rate ratio ~ 2.37, p = [fill in].
#6-azacytidine vs control: rate ratio ~ 0.92, p = [fill in].

emm_darting <- emmeans(component_models[["Darting"]], ~ Treatment, type = "response")
emm_darting
pairs(emm_darting, reverse = TRUE)
#Model-estimated darting was 0.30 events per observation in control boxes, 0.34 under
#6-azacytidine, and 0.71 under Decitabine. Decitabine vs control: rate ratio ~ 2.32,
#p = [fill in]. Decitabine vs 6-azacytidine: rate ratio ~ 2.10, p = [fill in].
#6-azacytidine vs control: rate ratio ~ 1.11, p = [fill in].

emm_humming <- emmeans(component_models[["Humming"]], ~ Treatment, type = "response")
emm_humming
pairs(emm_humming, reverse = TRUE)
#Model-estimated humming was 0.76 events per observation in control boxes, 1.23 under
#6-azacytidine, and 1.60 under Decitabine. Decitabine vs control: rate ratio ~ 2.11,
#p = [fill in]. Decitabine vs 6-azacytidine: rate ratio ~ 1.30, p = [fill in].
#6-azacytidine vs control: rate ratio ~ 1.62, p = [fill in] -- unlike attacking and
#darting, 6-azacytidine looks like it may differ from control here.

component_emm_df <- bind_rows(lapply(names(component_emm), function(nm) {
  as.data.frame(component_emm[[nm]]) %>%
    transmute(Behaviour = nm, Treatment, estimate = response,
              lower = asymp.LCL, upper = asymp.UCL)
}))

component_pairs_df <- bind_rows(lapply(names(component_pairs), function(nm) {
  as.data.frame(confint(component_pairs[[nm]])) %>%
    mutate(Behaviour = nm) %>%
    relocate(Behaviour)
}))

cat("\n--- Component marginal means (response scale) ---\n")
print(component_emm_df)

cat("\n--- Component pairwise contrasts (rate ratio scale, Tukey-adjusted) ---\n")
print(component_pairs_df)

## Read component_pairs_df alongside the Total pairs() output above:
## if all three components show the same Decitabine > 6-azacytidine >
## control ordering as Total, that's the case for treating them as
## biologically equivalent and the composite as a reasonable summary.
## If one component decouples from that pattern (e.g. Decitabine
## raises attacking specifically but not darting/humming, or vice
## versa), that's worth a sentence in the main text and the full
## component table belongs in the supplement.

## ---- Raw magnitude check: does one behaviour dominate Total? ----
## If one component contributes most of the counts, Total is really
## just a proxy for that one behaviour rather than a genuine composite.
component_totals <- data.frame(
  Behaviour   = component_vars,
  total_count = sapply(component_vars, function(v) sum(comp_data[[v]], na.rm = TRUE))
)
component_totals$proportion_of_sum <- component_totals$total_count / sum(component_totals$total_count)

cat("\n--- Raw totals across all observations, each component ---\n")
print(component_totals)


####################################
# --- Raw data dots + model lines/ribbon (population-level) ---

# Model predictions as you already built them
day_seq   <- sort(unique(comp_data$DayNum))
time_vals <- sort(unique(comp_data$Time))

newdat2 <- expand.grid(
  Treatment = levels(comp_data$Treatment),
  DayNum    = day_seq,
  Time      = time_vals
)

pred_link2 <- predict(
  m_total_spline,
  newdata = newdat2,
  type    = "link",
  se.fit  = TRUE,
  re.form = NA
)

newdat2 <- newdat2 %>%
  mutate(
    eta = pred_link2$fit,
    se  = pred_link2$se.fit,
    mu  = exp(eta),
    lcl = exp(eta - 1.96 * se),
    ucl = exp(eta + 1.96 * se)
  )

newdat_avg <- newdat2 %>%
  group_by(Treatment, DayNum) %>%
  summarise(
    mu  = mean(mu),
    lcl = mean(lcl),
    ucl = mean(ucl),
    .groups = "drop"
  )



# --- Marginal treatment means (right-hand inset) ---
emm_treat <- emmeans(m_total_spline, ~ Treatment, type = "response")
emm_df <- as.data.frame(emm_treat) %>%
  transmute(
    Treatment,
    mu  = response,
    lcl = asymp.LCL,
    ucl = asymp.UCL
  )

# --- x placement for the inset panel ---
x_min <- min(comp_data$DayNum, na.rm = TRUE)
x_max <- max(comp_data$DayNum, na.rm = TRUE)
x_rng <- x_max - x_min

x_inset_left  <- x_max + 0.03 * x_rng
x_inset_right <- x_max + 0.18 * x_rng
x_inset_mid   <- (x_inset_left + x_inset_right) / 2

emm_df$x <- x_inset_mid

# --- y helper for label placement ---
y_top <- max(c(comp_data$Total, newdat_avg$ucl, emm_df$ucl), na.rm = TRUE)

# --- colours (edit if you want different ones) ---
treat_cols <- c(
  "Control"       = "#333333",
  "6-azacytidine"  = "#1b9e77",
  "Decitabine"     = "#d95f02"
)

p <- ggplot() +

  # Inset background (white box)
  annotate(
    "rect",
    xmin = x_inset_left, xmax = x_inset_right,
    ymin = -Inf, ymax = Inf,
    fill = "white", colour = NA
  ) +

  # Separator line to mark inset boundary
  geom_vline(xintercept = x_inset_left, linewidth = 0.4, colour = "grey70") +

  # Model CI ribbon (behind)
  geom_ribbon(
    data = newdat_avg,
    aes(x = DayNum, ymin = lcl, ymax = ucl, fill = Treatment),
    alpha = 0.18,
    colour = NA
  ) +

  # Model mean line
  geom_line(
    data = newdat_avg,
    aes(x = DayNum, y = mu, colour = Treatment),
    linewidth = 1.0
  ) +

  # Raw observations
  geom_point(
    data = comp_data,
    aes(x = DayNum, y = Total, colour = Treatment),
    alpha = 0.18,
    size = 1.2,
    position = position_jitter(width = 0.10, height = 0)
  ) +

  # Marginal treatment means (95% CI) in inset
  geom_pointrange(
    data = emm_df,
    aes(x = x, y = mu, ymin = lcl, ymax = ucl, colour = Treatment),
    linewidth = 0.8
  ) +

  # Label for inset: left-aligned inside inset (won't overlap separator)
  annotate(
    "text",
    x = x_inset_left + 0.01 * x_rng,
    y = y_top * 1.04,
    label = "Marginal mean\n(95% CI)",
    hjust = 0,
    size = 3.2
  ) +

  scale_colour_manual(name = "Treatment", values = treat_cols) +
  scale_fill_manual(name = "Treatment", values = treat_cols) +

  labs(
    x = "Experimental day (within batch)",
    y = "Aggressive interactions (count)"
  ) +

  coord_cartesian(
    xlim = c(x_min, x_inset_right),
    ylim = c(0, NA),     # start y-axis at 0; remove if you don't want this
    clip = "off"
  ) +

  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.margin = margin(5.5, 10, 5.5, 5.5)
  )

p

ggsave(file.path(output_dir, "aggression.pdf"), p, width = 6.8, height = 4.6, units = "in")


## =========================================================
## NEW: same plot, one per component behaviour
## =========================================================
## Identical layout to aggression.pdf (raw points + spline model
## line/ribbon + marginal-mean inset), generalised to take a response
## column and its fitted component model so Attack/Darting/Humming can
## each get their own figure without copy-pasting the block three
## times. Reuses day_seq, time_vals and treat_cols already defined
## above.

build_component_plot <- function(response_var, model, ylab, out_file) {

  newdat <- expand.grid(
    Treatment = levels(comp_data$Treatment),
    DayNum    = day_seq,
    Time      = time_vals
  )

  pred_link <- predict(
    model,
    newdata = newdat,
    type    = "link",
    se.fit  = TRUE,
    re.form = NA
  )

  newdat <- newdat %>%
    mutate(
      eta = pred_link$fit,
      se  = pred_link$se.fit,
      mu  = exp(eta),
      lcl = exp(eta - 1.96 * se),
      ucl = exp(eta + 1.96 * se)
    )

  newdat_avg <- newdat %>%
    group_by(Treatment, DayNum) %>%
    summarise(mu = mean(mu), lcl = mean(lcl), ucl = mean(ucl), .groups = "drop")

  emm_resp <- emmeans(model, ~ Treatment, type = "response")
  emm_df_comp <- as.data.frame(emm_resp) %>%
    transmute(Treatment, mu = response, lcl = asymp.LCL, ucl = asymp.UCL)

  x_min <- min(comp_data$DayNum, na.rm = TRUE)
  x_max <- max(comp_data$DayNum, na.rm = TRUE)
  x_rng <- x_max - x_min

  x_inset_left  <- x_max + 0.03 * x_rng
  x_inset_right <- x_max + 0.18 * x_rng
  x_inset_mid   <- (x_inset_left + x_inset_right) / 2

  emm_df_comp$x <- x_inset_mid

  y_top <- max(c(comp_data[[response_var]], newdat_avg$ucl, emm_df_comp$ucl), na.rm = TRUE)

  p_comp <- ggplot() +

    annotate(
      "rect",
      xmin = x_inset_left, xmax = x_inset_right,
      ymin = -Inf, ymax = Inf,
      fill = "white", colour = NA
    ) +

    geom_vline(xintercept = x_inset_left, linewidth = 0.4, colour = "grey70") +

    geom_ribbon(
      data = newdat_avg,
      aes(x = DayNum, ymin = lcl, ymax = ucl, fill = Treatment),
      alpha = 0.18,
      colour = NA
    ) +

    geom_line(
      data = newdat_avg,
      aes(x = DayNum, y = mu, colour = Treatment),
      linewidth = 1.0
    ) +

    geom_point(
      data = comp_data,
      aes(x = DayNum, y = .data[[response_var]], colour = Treatment),
      alpha = 0.18,
      size = 1.2,
      position = position_jitter(width = 0.10, height = 0)
    ) +

    geom_pointrange(
      data = emm_df_comp,
      aes(x = x, y = mu, ymin = lcl, ymax = ucl, colour = Treatment),
      linewidth = 0.8
    ) +

    annotate(
      "text",
      x = x_inset_left + 0.01 * x_rng,
      y = y_top * 1.04,
      label = "Marginal mean\n(95% CI)",
      hjust = 0,
      size = 3.2
    ) +

    scale_colour_manual(name = "Treatment", values = treat_cols) +
    scale_fill_manual(name = "Treatment", values = treat_cols) +

    labs(
      x = "Experimental day (within batch)",
      y = ylab
    ) +

    coord_cartesian(
      xlim = c(x_min, x_inset_right),
      ylim = c(0, NA),
      clip = "off"
    ) +

    theme_classic(base_size = 12) +
    theme(
      legend.position = "bottom",
      plot.margin = margin(5.5, 10, 5.5, 5.5)
    )

  ggsave(file.path(output_dir, out_file), p_comp, width = 6.8, height = 4.6, units = "in")
  p_comp
}

p_attack <- build_component_plot(
  "Attack", component_models[["Attack"]],
  "Attacking behaviour (count)", "attacking.pdf"
)
p_attack

p_darting <- build_component_plot(
  "Darting", component_models[["Darting"]],
  "Darting behaviour (count)", "darting.pdf"
)
p_darting

p_humming <- build_component_plot(
  "Humming", component_models[["Humming"]],
  "Humming behaviour (count)", "humming.pdf"
)
p_humming


## =========================================================
## NEW: combine the three component plots into one supplementary figure
## =========================================================
library(patchwork)

supp_fig_components <- p_attack + p_darting + p_humming +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

supp_fig_components

ggsave(
  file.path(output_dir, "supp_fig_component_behaviours.pdf"),
  supp_fig_components,
  width  = 10,
  height = 4
)
