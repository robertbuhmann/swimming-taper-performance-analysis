# =============================================================================
# Swimming Taper Study
# Script 7: Secondary Aim — Change in Monitoring Variables Across the Taper
# Input:    taper_period_data.csv
# Output:   secondary_aim_results.txt, plots_lme_*.png
#
# Aim: Determine whether HRV (ln rMSSD), sleep quality, sleep duration,
# and mental fatigue differed between the pre-taper and taper periods,
# using linear mixed effects models with random intercepts for athlete.
# =============================================================================

library(tidyverse)
library(lme4)
library(lmerTest)   # provides p-values for lmer via Satterthwaite approximation
library(emmeans)    # estimated marginal means and contrasts
library(gridExtra)
library(lubridate)
library(robustlmm)

theme_set(
  theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey92"),
      legend.position  = "bottom"
    )
)


# =============================================================================
# 1. LOAD & PREPARE DATA
# =============================================================================

df <- read.csv("taper_period_analysis.csv") %>%
  mutate(
    date = as.Date(date, format = "%d/%m/%Y"),
    
    weekday = factor(
      weekday,
      levels = c("Monday", "Tuesday", "Wednesday",
                 "Thursday", "Friday", "Saturday", "Sunday")
    ),
    
    ln_rMSSD = log(rMSSD),
    sleep_quality_num = as.numeric(sleep_quality),
    
    sleep_time = if_else(sleep_time > 24 | sleep_time == 0,
                         NA_real_, sleep_time),
    sleep_quality_num = if_else(sleep_quality_num == 0,
                                NA_real_, sleep_quality_num),
    
    mental_fatigue_score = suppressWarnings(as.numeric(mental_fatigue_score)),
    
    week_of_year = isoweek(date),
    year = isoyear(date)
  )

weekly_summary <- df %>%
  group_by(week_of_year, athlete, taper, group, priority_event_1) %>%
  summarise(
    mean_ln_rMSSD = mean(ln_rMSSD, na.rm = TRUE),
    mean_sleep_time = mean(sleep_time, na.rm = TRUE),
    mean_sleep_qual = mean(sleep_quality, na.rm = TRUE),
    mean_mental_fatigue = mean(mental_fatigue_score, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

# Quick overview
cat("=== Data overview ===\n")
cat(sprintf("Total rows:  %d\n", nrow(weekly_summary)))
cat(sprintf("Athletes:    %d\n", n_distinct(weekly_summary$athlete)))
cat(sprintf("Date range:  %s to %s\n",
            min(df$date, na.rm = TRUE), max(weekly_summary$week_of_year, na.rm = TRUE)))
cat("\nPeriod breakdown:\n")
print(table(weekly_summary$taper, useNA = "ifany"))

cat("\nMissingness per outcome:\n")
weekly_summary %>%
  summarise(
    n_ln_rMSSD      = sum(!is.na(mean_ln_rMSSD)),
    n_sleep_quality  = sum(!is.na(mean_sleep_qual)),
    n_sleep_duration = sum(!is.na(mean_sleep_time)),
    n_mental_fatigue = sum(!is.na(mean_mental_fatigue)),
    miss_ln_rMSSD      = sum(is.na(mean_ln_rMSSD)),
    miss_sleep_quality  = sum(is.na(mean_sleep_qual)),
    miss_sleep_duration = sum(is.na(mean_sleep_time)),
    miss_mental_fatigue = sum(is.na(mean_mental_fatigue))
  ) %>%
  print()

# Check levels of all factor variables in the model
cat("--- period levels ---\n")
print(table(weekly_summary$taper, useNA = "ifany"))

cat("\n--- weekday levels ---\n")
print(table(weekly_summary$week_of_year, useNA = "ifany"))

cat("\n--- athlete levels ---\n")
cat(sprintf("%d unique athletes\n", n_distinct(weekly_summary$athlete)))

cat("\n--- taper column (raw) ---\n")
print(table(weekly_summary$taper, useNA = "ifany"))

cat("\n--- taper_flag column (raw) ---\n")
print(table(weekly_summary$taper, useNA = "ifany"))

# =============================================================================
# 2. FIT LINEAR MIXED EFFECTS MODELS
# =============================================================================
# One model per outcome:
#   outcome ~ period + weekday + (1 | athlete)
#
# Random intercept for athlete accounts for between-athlete baseline
# differences, making each athlete their own control.
# Weekday is included to account for within-week training load variation.
# Period (Pre-Taper vs Taper) is the predictor of interest.
# =============================================================================

fit_lme <- function(outcome_var, data) {
  formula <- as.formula(
    paste(outcome_var, "~ taper  + group*week_of_year + (1 | athlete)")
  )
  model <- tryCatch(
    lmer(formula, data = data, REML = TRUE),
    error = function(e) {
      message(sprintf("Model failed for %s: %s", outcome_var, e$message))
      NULL
    }
  )
  model
}

outcomes <- list(
  mean_ln_rMSSD             = "HRV - ln(rMSSD)",
  mean_sleep_qual    = "Sleep quality",
  mean_sleep_time           = "Sleep duration (hrs)",
  mean_mental_fatigue = "Mental fatigue score"
)

models <- map(names(outcomes), ~fit_lme(.x, weekly_summary))
names(models) <- names(outcomes)

# Check all models fitted
cat("\n=== Model convergence check ===\n")
walk2(models, names(models), function(m, nm) {
  cat(sprintf("%-25s %s\n", nm,
              if (is.null(m)) "FAILED" else "OK"))
})


# =============================================================================
# 3. EXTRACT RESULTS
# =============================================================================

# ── 3a. Random effects: Intraclass correlation (ICC) ─────────────────────────
# ICC = proportion of total variance explained by between-athlete differences
# Justifies the use of random intercepts

cat("\n=== ICC (between-athlete variance / total variance) ===\n")
icc_results <- map_df(names(models), function(nm) {
  m <- models[[nm]]
  if (is.null(m)) return(tibble(outcome = nm, icc = NA_real_))
  vc   <- as.data.frame(VarCorr(m))
  var_athlete   <- vc$vcov[vc$grp == "athlete"]
  var_residual  <- vc$vcov[vc$grp == "Residual"]
  icc <- var_athlete / (var_athlete + var_residual)
  tibble(
    outcome = outcomes[[nm]],
    var_athlete  = round(var_athlete, 4),
    var_residual = round(var_residual, 4),
    icc          = round(icc, 3)
  )
})
print(icc_results)
cat("ICC > 0.1 supports use of random intercept model\n")

# ── 3b. Fixed effects: period coefficient ────────────────────────────────────

cat("\n=== Fixed effects (period: Taper vs Pre-Taper) ===\n")
fixed_results <- map_df(names(models), function(nm) {
  m <- models[[nm]]
  if (is.null(m)) return(NULL)
  s    <- summary(m)
  coef_taper <- s$coefficients["taperY", , drop = FALSE]
  tibble(
    outcome    = outcomes[[nm]],
    beta       = round(coef_taper[, "Estimate"], 3),
    se         = round(coef_taper[, "Std. Error"], 3),
    df         = round(coef_taper[, "df"], 1),
    t          = round(coef_taper[, "t value"], 2),
    p          = round(coef_taper[, "Pr(>|t|)"], 3)
  )
})

print(fixed_results)

# ── 3c. Estimated marginal means and contrasts ────────────────────────────────
# emmeans gives model-adjusted means per period, averaging over weekday
# The contrast (Taper - Pre-Taper) is what goes in the results section

cat("\n=== Estimated marginal means and contrasts ===\n")
emm_results <- map_df(names(models), function(nm) {
  m <- models[[nm]]
  if (is.null(m)) return(NULL)
  
  emm  <- emmeans(m, ~ taper)
  cont <- contrast(emm, method = "revpairwise")  # Taper - Pre-Taper
  cont_df <- as.data.frame(cont)
  emm_df  <- as.data.frame(emm)
  
  tibble(
    outcome        = outcomes[[nm]],
    mean_pretaper  = round(emm_df$emmean[emm_df$taper == "N"], 3),
    mean_taper     = round(emm_df$emmean[emm_df$taper == "Y"], 3),
    difference     = round(cont_df$estimate, 3),      # Taper - Pre-Taper
    se_diff        = round(cont_df$SE, 3),
    ci_lower       = round(cont_df$estimate - 1.96 * cont_df$SE, 3),
    ci_upper       = round(cont_df$estimate + 1.96 * cont_df$SE, 3),
    t              = round(cont_df$t.ratio, 2),
    df             = round(cont_df$df, 1),
    p              = round(cont_df$p.value, 3),
    direction      = if_else(cont_df$estimate > 0, "increased", "decreased"),
    significant    = cont_df$p.value < 0.05
  )
})

cat("\n--- Estimated marginal means ---\n")
print(emm_results %>%
        select(outcome, mean_pretaper, mean_taper, difference,
               ci_lower, ci_upper, p, direction, significant))

# ── 3d. Full model summaries ──────────────────────────────────────────────────

cat("\n=== Full model summaries ===\n")
walk2(models, names(models), function(m, nm) {
  if (is.null(m)) return()
  cat(sprintf("\n--- %s ---\n", outcomes[[nm]]))
  print(summary(m))
})

# =============================================================================
# Secondary Aim: LME Model Diagnostics
# =============================================================================

library(lattice)   # qqmath for lmer residuals

# Named vector for clean plot titles
outcome_labels <- c(
  ln_rMSSD             = "HRV - ln(rMSSD)",
  sleep_quality_num    = "Sleep quality",
  sleep_time           = "Sleep duration (hrs)",
  mental_fatigue_score = "Mental fatigue score"
)

png("diag_lme_all.png", width = 1600, height = 1400, res = 130)
par(mfrow = c(4, 4),
    mar   = c(4, 4, 3, 1))

walk(names(models), function(nm) {
  
  m     <- models[[nm]]
  label <- outcome_labels[nm]
  
  if (is.null(m)) {
    plot.new(); title(paste(label, "— model failed"))
    return()
  }
  
  resids  <- residuals(m)
  fitted  <- fitted(m)
  hat_raw <- hatvalues(m)
  
  # ── Plot 1: Residuals vs Fitted (linearity + homoscedasticity) ─────────────
  plot(fitted, resids,
       xlab = "Fitted values",
       ylab = "Residuals",
       main = paste(label, "\nResid vs Fitted"),
       pch  = 19, col = "#4472C480", cex = 0.7)
  abline(h = 0, lty = 2, col = "red")
  lines(lowess(fitted, resids), col = "#ED7D31", lwd = 2)
  
  # ── Plot 2: Q-Q plot of residuals (normality) ──────────────────────────────
  qqnorm(resids,
         main = paste(label, "\nQ-Q Residuals"),
         pch  = 19, col = "#4472C480", cex = 0.7)
  qqline(resids, col = "red", lwd = 2)
  
  # ── Plot 3: Scale-Location (homoscedasticity) ──────────────────────────────
  plot(fitted, sqrt(abs(resids)),
       xlab = "Fitted values",
       ylab = expression(sqrt("|Residuals|")),
       main = paste(label, "\nScale-Location"),
       pch  = 19, col = "#4472C480", cex = 0.7)
  lines(lowess(fitted, sqrt(abs(resids))), col = "#ED7D31", lwd = 2)
  
  # ── Plot 4: Leverage / influential observations ────────────────────────────
  # For LME, use standardised residuals vs leverage (hat values)
  # Flag points exceeding the median leverage + 3*IQR
  std_res      <- resids / sd(resids)
  lev_threshold <- median(hat_raw) + 3 * IQR(hat_raw)
  flagged      <- hat_raw > lev_threshold | abs(std_res) > 2
  
  plot(hat_raw, std_res,
       xlab = "Leverage (hat value)",
       ylab = "Standardised residual",
       main = paste(label, "\nLeverage"),
       pch  = 19,
       col  = ifelse(flagged, "#E63946", "#4472C480"),
       cex  = 0.7)
  abline(h = c(-2, 2), lty = 2, col = "grey50")
  abline(v = lev_threshold, lty = 2, col = "#ED7D31")
  
  # Label flagged points with observation number
  if (any(flagged)) {
    text(hat_raw[flagged], std_res[flagged],
         labels = which(flagged),
         cex = 0.6, pos = 3, col = "#E63946")
  }
})

par(mfrow = c(1, 1))
dev.off()
cat("Saved: diag_lme_all.png\n")


# ── Shapiro-Wilk on residuals for each model ──────────────────────────────────

cat("\n--- Shapiro-Wilk test on residuals ---\n")
map_df(names(models), function(nm) {
  m <- models[[nm]]
  if (is.null(m)) return(tibble(outcome = nm, W = NA, p = NA))
  sw <- shapiro.test(residuals(m)[seq_len(min(length(residuals(m)), 5000))])
  tibble(
    outcome = outcome_labels[nm],
    W       = round(sw$statistic, 4),
    p       = round(sw$p.value, 4),
    normal  = if_else(sw$p.value > 0.05, "Yes", "No")
  )
}) %>% print()


# ── Q-Q plot of random effects (normality of random intercepts) ───────────────

png("diag_lme_random_effects.png", width = 1200, height = 400, res = 130)
par(mfrow = c(1, 4), mar = c(4, 4, 3, 1))

walk(names(models), function(nm) {
  m <- models[[nm]]
  if (is.null(m)) { plot.new(); return() }
  re <- ranef(m)$athlete[, 1]
  qqnorm(re,
         main = paste(outcome_labels[nm], "\nRandom intercepts Q-Q"),
         pch  = 19, col = "#4472C4", cex = 0.8)
  qqline(re, col = "red", lwd = 2)
})

par(mfrow = c(1, 1))
dev.off()
cat("Saved: diag_lme_random_effects.png\n")

# =============================================================================
# 4. GENERATE RESULTS TEXT FOR MANUSCRIPT
# =============================================================================

n_athletes   <- n_distinct(df$athlete)
n_obs        <- nrow(df)

generate_result_sentence <- function(row) {
  direction_word <- if_else(row$significant,
                            row$direction,
                            paste("did not significantly change")
  )
  
  sig_note <- if_else(row$significant, "", " (non-significant)")
  
  sprintf(
    "%s %s during the taper period relative to the pre-taper period (estimated marginal mean difference: %.3f, 95%% CI: %.3f, %.3f, p = %.3f)%s.",
    row$outcome,
    direction_word,
    row$difference,
    row$ci_lower,
    row$ci_upper,
    row$p,
    sig_note
  )
}

cat("\n\n")
cat("=================================================================\n")
cat("MANUSCRIPT TEXT — Secondary Aim Results\n")
cat("=================================================================\n\n")

cat(sprintf(
  "Linear mixed effects models with random intercepts for athlete were used
to examine whether HRV (ln rMSSD), sleep quality, sleep duration, and
mental fatigue differed between the pre-taper and taper periods across the
full daily-level dataset (n = %d athletes, %d total observations). Weekday
was included as a covariate to account for within-week variation in
training load.\n\n",
  n_athletes, n_obs
))

walk(1:nrow(emm_results), function(i) {
  cat(generate_result_sentence(emm_results[i, ]), "\n\n")
})

cat("=================================================================\n\n")


# =============================================================================
# 5. VISUALISATION
# =============================================================================

# ── 5a. Estimated marginal means plot ────────────────────────────────────────

emm_plot_data <- map_df(names(models), function(nm) {
  m <- models[[nm]]
  if (is.null(m)) return(NULL)
  emm <- as.data.frame(emmeans(m, ~ taper))
  emm$outcome <- outcomes[[nm]]
  emm
})

p_emm <- ggplot(emm_plot_data,
                aes(x = taper, y = emmean, colour = taper, group = outcome)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                width = 0.2, linewidth = 0.8) +
  geom_line(aes(group = 1), colour = "grey50", linewidth = 0.5,
            linetype = "dashed") +
  scale_colour_manual(values = c("Pre-Taper" = "#4472C4",
                                 "Taper"     = "#ED7D31")) +
  facet_wrap(~outcome, scales = "free_y", ncol = 2) +
  labs(
    title    = "Estimated Marginal Means: Pre-Taper vs Taper",
    subtitle = "Adjusted for weekday and between-athlete variation",
    x        = NULL,
    y        = "Estimated marginal mean (95% CI)",
    colour   = "Period"
  ) +
  theme(legend.position = "bottom")

ggsave("plots_lme_emm.png", p_emm, width = 10, height = 8, dpi = 300)
cat("Saved: plots_lme_emm.png\n")

# ── 5b. Individual trajectory plots (spaghetti) ───────────────────────────────

make_spaghetti <- function(var, label) {
  ggplot(df %>% filter(!is.na(.data[[var]])),
         aes(x = date, y = .data[[var]],
             group = athlete, colour = athlete)) +
    geom_line(alpha = 0.4, linewidth = 0.5) +
    geom_smooth(aes(group = 1), method = "loess", span = 0.4,
                colour = "black", linewidth = 1.2, se = FALSE) +
    facet_wrap(~period, scales = "free_x") +
    labs(title = paste("Individual trajectories:", label),
         x = "Date", y = label) +
    theme(legend.position = "none")
}

p_hrv_sp   <- make_spaghetti("ln_rMSSD",          "HRV - ln(rMSSD)")
p_sleep_sp <- make_spaghetti("sleep_quality_num",  "Sleep quality")
p_dur_sp   <- make_spaghetti("sleep_time",         "Sleep duration (hrs)")
p_mf_sp    <- make_spaghetti("mental_fatigue_score","Mental fatigue score")

png("plots_lme_spaghetti.png", width = 1200, height = 1400, res = 120)
gridExtra::grid.arrange(p_hrv_sp, p_sleep_sp, p_dur_sp, p_mf_sp, ncol = 1)
dev.off()
cat("Saved: plots_lme_spaghetti.png\n")

# ── 5c. Box plots with individual data points ─────────────────────────────────

bp_long <- df %>%
  select(athlete, period, ln_rMSSD, sleep_quality_num,
         sleep_time, mental_fatigue_score) %>%
  pivot_longer(
    cols      = c(ln_rMSSD, sleep_quality_num,
                  sleep_time, mental_fatigue_score),
    names_to  = "variable",
    values_to = "value"
  ) %>%
  filter(!is.na(value)) %>%
  mutate(variable = case_match(variable,
                               "ln_rMSSD"             ~ "HRV - ln(rMSSD)",
                               "sleep_quality_num"    ~ "Sleep quality",
                               "sleep_time"           ~ "Sleep duration (hrs)",
                               "mental_fatigue_score" ~ "Mental fatigue score",
                               .default = variable
  ))

# Add p-value annotations from emm_results
p_annotations <- emm_results %>%
  mutate(
    variable = outcome,
    label    = paste0("p = ", formatC(p, digits = 3, format = "f"))
  )

p_boxplots <- ggplot(bp_long,
                     aes(x = period, y = value, fill = period)) +
  geom_boxplot(alpha = 0.6, outlier.shape = 21, width = 0.45) +
  geom_jitter(aes(colour = period), width = 0.12,
              alpha = 0.25, size = 0.9) +
  geom_text(
    data = p_annotations %>%
      left_join(bp_long %>%
                  group_by(variable) %>%
                  summarise(y_pos = max(value, na.rm = TRUE) * 1.05,
                            .groups = "drop"),
                by = "variable"),
    aes(x = 1.5, y = y_pos, label = label),
    inherit.aes = FALSE, size = 3.2
  ) +
  scale_fill_manual(values  = c("Pre-Taper" = "#4472C4",
                                "Taper"     = "#ED7D31")) +
  scale_colour_manual(values = c("Pre-Taper" = "#4472C4",
                                 "Taper"     = "#ED7D31")) +
  facet_wrap(~variable, scales = "free_y", ncol = 2) +
  labs(
    title = "Monitoring Variables: Pre-Taper vs Taper",
    x     = NULL,
    y     = "Value",
    fill  = "Period"
  ) +
  theme(legend.position = "bottom",
        legend.title    = element_blank()) +
  guides(colour = "none")

ggsave("plots_lme_boxplots.png", p_boxplots,
       width = 10, height = 8, dpi = 300)
cat("Saved: plots_lme_boxplots.png\n")


# =============================================================================
# 6. SAVE RESULTS
# =============================================================================

write.csv(emm_results,   "secondary_aim_emm_results.csv",   row.names = FALSE)
write.csv(fixed_results, "secondary_aim_fixed_effects.csv",  row.names = FALSE)
write.csv(icc_results,   "secondary_aim_icc.csv",            row.names = FALSE)
cat("Saved: secondary_aim_emm_results.csv\n")
cat("Saved: secondary_aim_fixed_effects.csv\n")
cat("Saved: secondary_aim_icc.csv\n")

# Full text report
sink("secondary_aim_results.txt")

cat("=== Secondary Aim: LME Results ===\n\n")
cat(sprintf("n athletes:     %d\n", n_athletes))
cat(sprintf("n observations: %d\n\n", n_obs))

cat("--- ICC (between-athlete variance) ---\n")
print(icc_results)

cat("\n--- Estimated marginal means and contrasts (Taper - Pre-Taper) ---\n")
print(as.data.frame(emm_results))

cat("\n--- Manuscript text ---\n\n")
cat(sprintf(
  "Linear mixed effects models with random intercepts for athlete were used
to examine whether HRV (ln rMSSD), sleep quality, sleep duration, and
mental fatigue differed between the pre-taper and taper periods across the
full daily-level dataset (n = %d athletes, %d total observations). Weekday
was included as a covariate to account for within-week variation in
training load.\n\n",
  n_athletes, n_obs
))

walk(1:nrow(emm_results), function(i) {
  cat(generate_result_sentence(emm_results[i, ]), "\n\n")
})

sink()
cat("Saved: secondary_aim_results.txt\n")
cat("\n=== SCRIPT 7 COMPLETE ===\n")
