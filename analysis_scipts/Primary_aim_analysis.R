# ============================================================
# Swimming taper study
# Athlete-level modelling script with robust standard errors
# ============================================================

library(tidyverse)
library(car)
library(lmtest)
library(sandwich)

theme_set(
  theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey92"),
      legend.position = "bottom"
    )
)

# ============================================================
# 1. Load and prepare data
# ============================================================

df <- readRDS("swimming_training_load_imputed_2026-03-31.rds")

df <- df %>%
  mutate(
    date = as.Date(date),
    athlete = factor(athlete),
    group = factor(group),
    period = factor(
      if_else(taper == 1 | taper_flag == 1, "Taper", "Pre-Taper"),
      levels = c("Pre-Taper", "Taper")
    ),
    mental_fatigue_score = as.numeric(mental_fatigue_score),
    physical_fatigue_score = as.numeric(physical_fatigue_score),
    ln_rMSSD = log(rMSSD),
    sleep_quality_num = as.numeric(sleep_quality),
    perf_change = perf_change1,
    perf_change_s = perf_change1_s,
    perf_change_pct = per_change1_pct
  )

# ============================================================
# 2. Clean key variables
# ============================================================

df <- df %>%
  mutate(
    sleep_time = if_else(sleep_time > 24, NA_real_, sleep_time),
    sleep_quality_num = if_else(sleep_quality_num == 0, NA_real_, sleep_quality_num)
  )

cat("Rows:", nrow(df), "
")
cat("Athletes:", n_distinct(df$athlete), "
")
cat("Pre-Taper rows:", sum(df$period == "Pre-Taper", na.rm = TRUE), "
")
cat("Taper rows:", sum(df$period == "Taper", na.rm = TRUE), "
")

# ============================================================
# 3. Create athlete-level dataset
# ============================================================

df_athlete <- df %>%
  group_by(athlete) %>%
  summarise(
    taper_ln_rMSSD = mean(ln_rMSSD[period == "Taper"], na.rm = TRUE),
    taper_sleep_quality = mean(sleep_quality_num[period == "Taper"], na.rm = TRUE),
    taper_sleep_duration = mean(sleep_time[period == "Taper"], na.rm = TRUE),
    taper_mental_fatigue = mean(mental_fatigue_score[period == "Taper"], na.rm = TRUE),
    delta_ln_rMSSD =
      mean(ln_rMSSD[period == "Taper"], na.rm = TRUE) -
      mean(ln_rMSSD[period == "Pre-Taper"], na.rm = TRUE),
    delta_sleep_quality =
      mean(sleep_quality_num[period == "Taper"], na.rm = TRUE) -
      mean(sleep_quality_num[period == "Pre-Taper"], na.rm = TRUE),
    delta_mental_fatigue =
      mean(mental_fatigue_score[period == "Taper"], na.rm = TRUE) -
      mean(mental_fatigue_score[period == "Pre-Taper"], na.rm = TRUE),
    perf_change = first(perf_change),
    perf_change_s = first(perf_change_s),
    perf_change_pct = first(perf_change_pct),
    group = first(group),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.), NA, .)))

# ============================================================
# 4. Scale taper predictors for the primary model
# ============================================================

df_athlete <- df_athlete %>%
  mutate(
    taper_ln_rMSSD_scaled = as.numeric(scale(taper_ln_rMSSD)),
    taper_sleep_quality_scaled = as.numeric(scale(taper_sleep_quality)),
    taper_sleep_duration_scaled = as.numeric(scale(taper_sleep_duration)),
    taper_mental_fatigue_scaled = as.numeric(scale(taper_mental_fatigue))
  )

cat("Athletes in athlete-level dataset:", nrow(df_athlete), "
")

# ============================================================
# 5. Build model datasets and fit models
# ============================================================

m_primary_data <- df_athlete %>%
  select(
    athlete,
    perf_change_pct,
    taper_ln_rMSSD_scaled,
    taper_sleep_quality_scaled,
    taper_sleep_duration_scaled,
    taper_mental_fatigue_scaled
  )

m_delta_data <- df_athlete %>%
  select(
    athlete,
    perf_change_pct,
    delta_ln_rMSSD,
    delta_sleep_quality,
    delta_mental_fatigue
  )

m_primary_data_used <- m_primary_data %>%
  filter(complete.cases(.))

m_delta_data_used <- m_delta_data %>%
  filter(complete.cases(.))

cat("Athletes in primary model:", nrow(m_primary_data_used), "
")
cat("Athletes in delta model:", nrow(m_delta_data_used), "
")

m_primary <- lm(
  perf_change_pct ~ taper_ln_rMSSD_scaled +
    taper_sleep_quality_scaled +
    taper_sleep_duration_scaled +
    taper_mental_fatigue_scaled,
  data = m_primary_data_used
)

m_delta <- lm(
  perf_change_pct ~ delta_ln_rMSSD +
    delta_sleep_quality +
    delta_mental_fatigue,
  data = m_delta_data_used
)

# ============================================================
# 6. Robust standard error helper (HC3)
# ============================================================

get_robust_results <- function(model, model_name) {
  vcov_hc3 <- sandwich::vcovHC(model, type = "HC3")
  robust_test <- lmtest::coeftest(model, vcov. = vcov_hc3)
  robust_ci <- lmtest::coefci(model, vcov. = vcov_hc3)

  robust_table <- data.frame(
    model = model_name,
    term = rownames(robust_test),
    estimate = robust_test[, "Estimate"],
    robust_se = robust_test[, "Std. Error"],
    statistic = robust_test[, "t value"],
    p_value = robust_test[, "Pr(>|t|)"],
    ci_lower = robust_ci[, 1],
    ci_upper = robust_ci[, 2],
    row.names = NULL
  )

  return(robust_table)
}

primary_robust_results <- get_robust_results(m_primary, "Primary model")
delta_robust_results <- get_robust_results(m_delta, "Delta model")

cat("
================ PRIMARY MODEL (ROBUST SE, HC3) ================
")
print(primary_robust_results)
cat("
Model fit:
")
print(summary(m_primary))

cat("
================ DELTA MODEL (ROBUST SE, HC3) ================
")
print(delta_robust_results)
cat("
Model fit:
")
print(summary(m_delta))

# ============================================================
# 7. Save outputs for diagnostics script
# ============================================================

saveRDS(m_primary, "m_primary_results.RDS")
saveRDS(m_delta, "m_delta_results.RDS")
saveRDS(m_primary_data_used, "m_primary_data_used.RDS")
saveRDS(m_delta_data_used, "m_delta_data_used.RDS")
write.csv(df_athlete, "df_athlete_model_data.csv", row.names = FALSE)
write.csv(m_primary_data_used, "m_primary_data_used.csv", row.names = FALSE)
write.csv(m_delta_data_used, "m_delta_data_used.csv", row.names = FALSE)
write.csv(primary_robust_results, "m_primary_robust_results.csv", row.names = FALSE)
write.csv(delta_robust_results, "m_delta_robust_results.csv", row.names = FALSE)

cat("
Saved files:
")
cat("- m_primary_results.RDS
")
cat("- m_delta_results.RDS
")
cat("- m_primary_data_used.RDS
")
cat("- m_delta_data_used.RDS
")
cat("- df_athlete_model_data.csv
")
cat("- m_primary_data_used.csv
")
cat("- m_delta_data_used.csv
")
cat("- m_primary_robust_results.csv
")
cat("- m_delta_robust_results.csv
")

