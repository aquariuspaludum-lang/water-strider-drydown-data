# ============================================================
# Stage-specific effects of habitat desiccation on larval
# developments in four Canadian water strider species
# Kishi & Spence (in preparation)
#
# Statistical analyses for manuscript tables and figures
# ============================================================
#
# Required packages:
#   install.packages(c("lme4","lmerTest","emmeans","multcomp","tidyverse"))
#
# Working directory should contain all CSV data files.
# ============================================================

library(lme4)
library(lmerTest)
library(emmeans)
library(multcomp)
library(tidyverse)

# ── Data loading ─────────────────────────────────────────────

np_files <- list(
  G_pingreensis = "Nymphal_period_G_pingreensis.csv",
  G_buenoi      = "Nymphal_period_G_buenoi.csv",
  A_remigis     = "Nymphal_period_A_remigis.csv",
  L_dissortis   = "Nymphal_period_L_dissortis.csv"
)

df_wing  <- read.csv("Wing_types.csv")
df_fm_bu <- read.csv("flight_muscle_G_buenoi.csv")
df_fm_re <- read.csv("flight_muscle_A_remigis.csv")
df_fm_di <- read.csv("flight_muscle_L_dissortis.csv")

preov_bu_lw <- read.csv("preoviposition_period_G_buenoi_long_wing.csv")
preov_bu_sw <- read.csv("preoviposition_period_G_buenoi_short_wing.csv")
preov_re_lw <- read.csv("preoviposition_period_A_remigis_long_wing.csv")
preov_re_sw <- read.csv("preoviposition_period_A_remigis_short_wing.csv")

# Starting N per condition (n per replicate case × number of replicates)
start_n <- list(
  G_pingreensis = 325,  # 65/case × 5 replicates
  G_buenoi      = 236,  # 59/case × 4 replicates
  A_remigis     = 207,  # 23/case × 9 replicates
  L_dissortis   = 165   # 33/case × 5 replicates
)


# ============================================================
# Analysis A: Total nymphal period
# GLMM: Gamma distribution, Replicate as random effect
# ============================================================

cat("\n====================================================\n")
cat("Analysis A: Nymphal period (GLMM)\n")
cat("====================================================\n")

results_np <- list()

for (sp in names(np_files)) {
  cat("\n---", sp, "---\n")

  df <- read.csv(np_files[[sp]]) %>%
    filter(Stage == "Adult") %>%
    mutate(Condition = factor(Condition),
           Replicate = factor(Replicate))

  model_full <- glmer(
    Days_to_stage ~ Condition + (1 | Replicate),
    family  = Gamma(link = "log"),
    data    = df,
    control = glmerControl(optimizer = "bobyqa",
                           optCtrl   = list(maxfun = 2e5))
  )

  model_null <- glmer(
    Days_to_stage ~ 1 + (1 | Replicate),
    family  = Gamma(link = "log"),
    data    = df,
    control = glmerControl(optimizer = "bobyqa")
  )

  cat("LRT (Condition effect):\n")
  print(anova(model_null, model_full))

  emm <- emmeans(model_full, ~ Condition)
  cat("\nPost-hoc vs Control (Dunnett):\n")
  print(summary(contrast(emm, method = "trt.vs.ctrl", ref = 1),
                adjust = "dunnett"))

  means <- df %>%
    group_by(Condition) %>%
    summarise(n = n(), mean = mean(Days_to_stage),
              se = sd(Days_to_stage) / sqrt(n()), .groups = "drop")
  cat("\nDescriptive stats (Table 2):\n")
  print(means)

  results_np[[sp]] <- list(model = model_full, emm = emm, means = means)
}


# ============================================================
# Analysis B: Stage-specific mortality
# GLM: Binomial (aggregate count data)
#
# Mortality data computed directly from nymphal period CSVs.
# Stage column in CSV is character: "2","3","4","5","Adult".
# Starting N per condition = n_per_case × n_replicates (see start_n).
# ============================================================

cat("\n====================================================\n")
cat("Analysis B: Stage-specific mortality (GLM, binomial)\n")
cat("====================================================\n")

species_codes <- c(
  G_pingreensis = "G. pingreensis",
  G_buenoi      = "G. buenoi",
  A_remigis     = "A. remigis",
  L_dissortis   = "L. dissortis"
)

results_mor <- list()

for (sp in names(np_files)) {
  cat("\n---", species_codes[sp], "---\n")

  df_np   <- read.csv(np_files[[sp]])
  N_start <- start_n[[sp]]

  # Count individuals reaching each stage per condition
  # Stage is character: "2","3","4","5","Adult"
  cnt <- function(stage_val) {
    tabulate(df_np$Condition[df_np$Stage == stage_val], nbins = 8)
  }
  n2     <- cnt("2")
  n3     <- cnt("3")
  n4     <- cnt("4")
  n5     <- cnt("5")
  nadult <- cnt("Adult")

  # Build mortality data frame (40 rows: 8 conditions × 5 instars)
  df_sp <- data.frame(
    Condition = factor(rep(1:8, 5)),
    Instar    = factor(rep(1:5, each = 8)),
    Survived  = c(n2,           n3,      n4,      n5,      nadult),
    Died      = c(N_start - n2, n2 - n3, n3 - n4, n4 - n5, n5 - nadult)
  )

  cat(sprintf("  N_start=%d, Adults (Cond1)=%d, Adults (Cond8)=%d\n",
              N_start, nadult[1], nadult[8]))

  # B-1: Overall additive model (Condition + Instar)
  cat("\n[B-1] Overall LRT (Condition + Instar):\n")
  model_add    <- glm(cbind(Died, Survived) ~ Condition + Instar,
                      family = binomial, data = df_sp)
  model_instar <- glm(cbind(Died, Survived) ~ Instar,
                      family = binomial, data = df_sp)
  model_null   <- glm(cbind(Died, Survived) ~ 1,
                      family = binomial, data = df_sp)

  cat("  Instar effect:\n")
  print(anova(model_null, model_instar, test = "LRT"))
  cat("  Condition effect (controlling for Instar):\n")
  print(anova(model_instar, model_add, test = "LRT"))

  # B-2: Per-instar Condition effect
  cat("\n[B-2] Per-instar Condition effect:\n")
  instar_results <- list()

  for (inst in 1:5) {
    df_i <- droplevels(df_sp[df_sp$Instar == inst, ])

    if (sum(df_i$Died) == 0) {
      cat(sprintf("  Instar %d: no deaths (skipped)\n", inst))
      next
    }

    m_full <- glm(cbind(Died, Survived) ~ Condition,
                  family = binomial, data = df_i)
    m_null <- glm(cbind(Died, Survived) ~ 1,
                  family = binomial, data = df_i)

    lrt  <- anova(m_null, m_full, test = "LRT")
    chi2 <- lrt$Deviance[2]
    pval <- lrt$`Pr(>Chi)`[2]

    cat(sprintf("  Instar %d: chi2=%.2f, df=%d, p=%.4f%s\n",
                inst, chi2, lrt$Df[2], pval,
                ifelse(pval < 0.001, " ***",
                       ifelse(pval < 0.01,  " **",
                              ifelse(pval < 0.05,  " *",  "")))))
    instar_results[[as.character(inst)]] <- m_full
  }

  # B-3: Post-hoc for significant instars
  if (length(instar_results) > 0) {
    cat("\n[B-3] Post-hoc vs Control (Dunnett):\n")
    for (inst_key in names(instar_results)) {
      emm    <- emmeans(instar_results[[inst_key]], ~ Condition)
      ph_sum <- summary(contrast(emm, method = "trt.vs.ctrl", ref = 1),
                        adjust = "dunnett")
      if (any(ph_sum$p.value < 0.05, na.rm = TRUE)) {
        cat(sprintf("\n  Instar %s:\n", inst_key))
        print(ph_sum)
      }
    }
  }

  # B-4: Mortality rate table (Table 3)
  df_sp$mort_rate <- round(df_sp$Died / (df_sp$Died + df_sp$Survived) * 100, 1)
  summary_tab <- reshape(df_sp[, c("Condition","Instar","mort_rate")],
                         idvar = "Condition", timevar = "Instar",
                         direction = "wide")
  names(summary_tab)[-1] <- paste0("Instar_", 1:5)
  cat("\nMortality rate (%) — Table 3:\n")
  print(summary_tab)

  results_mor[[sp]] <- list(model_add = model_add, instar_models = instar_results)
}


# ============================================================
# Analysis C: Wing morph proportion
# Logistic regression: Group × Sex
# (G. buenoi and A. remigis only)
# ============================================================

cat("\n====================================================\n")
cat("Analysis C: Wing morph (GLM, binomial)\n")
cat("====================================================\n")

df_wing_wide <- df_wing %>%
  mutate(LW = ifelse(Wing_type == "long_wing",  N_of_individuals, 0),
         SW = ifelse(Wing_type == "short_wing", N_of_individuals, 0)) %>%
  group_by(Species, Sex, Group) %>%
  summarise(LW = sum(LW), SW = sum(SW), .groups = "drop") %>%
  mutate(Group = factor(Group))

for (sp in c("G_buenoi", "A_remigis")) {
  cat("\n---", sp, "---\n")
  df_sp <- df_wing_wide %>% filter(Species == sp)

  model_full <- glm(cbind(LW, SW) ~ Group * Sex, family = binomial, data = df_sp)
  model_add  <- glm(cbind(LW, SW) ~ Group + Sex, family = binomial, data = df_sp)
  model_sex  <- glm(cbind(LW, SW) ~ Sex,         family = binomial, data = df_sp)
  model_grp  <- glm(cbind(LW, SW) ~ Group,       family = binomial, data = df_sp)

  cat("LRT — Group × Sex interaction:\n")
  print(anova(model_add, model_full, test = "LRT"))
  cat("LRT — Group effect:\n")
  print(anova(model_sex, model_add, test = "LRT"))
  cat("LRT — Sex effect:\n")
  print(anova(model_grp, model_add, test = "LRT"))

  emm <- emmeans(model_add, ~ Group)
  cat("\nGroup pairwise contrasts (Holm):\n")
  print(summary(contrast(emm, method = "pairwise"), adjust = "holm"))

  cat("\nLW proportion (%) — Table 4:\n")
  print(df_sp %>% mutate(LW_pct = round(LW / (LW + SW) * 100, 1)))
}


# ============================================================
# Analysis D: Preoviposition period
# GLM: Gamma distribution, Condition + Wing
# (G. buenoi and A. remigis only)
# ============================================================

cat("\n====================================================\n")
cat("Analysis D: Preoviposition period (GLM, Gamma)\n")
cat("====================================================\n")

df_preov_bu <- bind_rows(
  preov_bu_lw %>% mutate(Wing = "Long"),
  preov_bu_sw %>% mutate(Wing = "Short")
) %>% mutate(Condition = factor(Condition), Wing = factor(Wing))

df_preov_re <- bind_rows(
  preov_re_lw %>% mutate(Wing = "Long"),
  preov_re_sw %>% mutate(Wing = "Short")
) %>% mutate(Condition = factor(Condition), Wing = factor(Wing))

for (info in list(list(sp = "G_buenoi",  df = df_preov_bu),
                  list(sp = "A_remigis", df = df_preov_re))) {

  cat("\n---", info$sp, "---\n")
  df <- info$df

  model_full <- glm(Preoviposition ~ Condition * Wing,
                    family = Gamma(link = "log"), data = df)
  model_add  <- glm(Preoviposition ~ Condition + Wing,
                    family = Gamma(link = "log"), data = df)
  model_wing <- glm(Preoviposition ~ Wing,
                    family = Gamma(link = "log"), data = df)
  model_cond <- glm(Preoviposition ~ Condition,
                    family = Gamma(link = "log"), data = df)

  cat("LRT — Condition × Wing interaction:\n")
  print(anova(model_add, model_full, test = "LRT"))
  cat("LRT — Condition effect:\n")
  print(anova(model_wing, model_add, test = "LRT"))
  cat("LRT — Wing effect:\n")
  print(anova(model_cond, model_add, test = "LRT"))

  means <- df %>%
    group_by(Condition, Wing) %>%
    summarise(n = n(), mean = mean(Preoviposition),
              se = sd(Preoviposition) / sqrt(n()), .groups = "drop")
  cat("\nDescriptive stats (Table 5):\n")
  print(means)

  emm <- emmeans(model_add, ~ Condition)
  cat("\nCondition pairwise vs Control (Dunnett):\n")
  print(summary(contrast(emm, method = "trt.vs.ctrl", ref = 1),
                adjust = "dunnett"))
}


# ============================================================
# Analysis E: Flight muscle condition
# Chi-square test (condition effect on muscle category)
# ============================================================

cat("\n====================================================\n")
cat("Analysis E: Flight muscle (Chi-square)\n")
cat("====================================================\n")

fm_data <- list(G_buenoi = df_fm_bu, A_remigis = df_fm_re, L_dissortis = df_fm_di)

for (sp in names(fm_data)) {
  cat("\n---", sp, "---\n")
  df <- fm_data[[sp]]

  table_mat   <- as.matrix(df[, c("Developed","Intermediate","Histolysed")])
  rownames(table_mat) <- paste0("Cond", df$Condition)
  table_clean <- table_mat[, colSums(table_mat) > 0]  # drop all-zero columns

  cat("Contingency table:\n")
  print(table_clean)

  test <- chisq.test(table_clean)
  cat(sprintf("\nChi-square: X2=%.2f, df=%d, p=%.4f\n",
              test$statistic, test$parameter, test$p.value))

  cat("Proportions (%):\n")
  print(round(prop.table(table_clean, margin = 1) * 100, 1))
}


# ============================================================
# Table 2: Nymphal period summary (all species)
# ============================================================

cat("\n====================================================\n")
cat("Table 2: Nymphal period summary\n")
cat("====================================================\n")

table2 <- map_dfr(names(np_files), function(sp) {
  read.csv(np_files[[sp]]) %>%
    filter(Stage == "Adult") %>%
    group_by(Condition) %>%
    summarise(n    = n(),
              Mean = round(mean(Days_to_stage), 1),
              SE   = round(sd(Days_to_stage) / sqrt(n()), 1),
              .groups = "drop") %>%
    mutate(Species = sp)
}) %>% select(Species, Condition, n, Mean, SE)

print(table2, n = 32)

cat("\nCondition codes:\n")
cat("  1 = Control (water throughout)\n")
cat("  2 = Damp paper during 1st instar\n")
cat("  3 = Damp paper during 2nd instar\n")
cat("  4 = Damp paper during 3rd instar\n")
cat("  5 = Damp paper during 4th instar\n")
cat("  6 = Damp paper during early 5th instar\n")
cat("  7 = Damp paper during late 5th instar\n")
cat("  8 = Damp paper throughout (all instars)\n")


# ============================================================
# Analysis F: Flight muscle — inter-species comparison
# Chi-square / Fisher's exact:
# Maintained (Developed + Intermediate) vs Histolysed
# across G. buenoi, A. remigis, and L. dissortis
# ============================================================

cat("\n====================================================\n")
cat("Analysis F: Flight muscle — inter-species comparison\n")
cat("  (Maintained [Dev+Int] vs Histolysed, pooled across conditions)\n")
cat("====================================================\n")

# Pool all conditions per species: Maintained vs Histolysed
pool_maint <- function(df) {
  c(Maintained = sum(df$Developed) + sum(df$Intermediate),
    Histolysed  = sum(df$Histolysed))
}

tab_intersp <- rbind(
  G_buenoi    = pool_maint(df_fm_bu),
  A_remigis   = pool_maint(df_fm_re),
  L_dissortis = pool_maint(df_fm_di)
)

cat("\nContingency table (Maintained vs Histolysed, pooled across conditions):\n")
print(tab_intersp)

cat("\nHistolysed proportion per species (%):\n")
print(round(tab_intersp[, "Histolysed"] / rowSums(tab_intersp) * 100, 1))

cat("\nMaintained proportion per species (%):\n")
print(round(tab_intersp[, "Maintained"] / rowSums(tab_intersp) * 100, 1))

# Overall 3-species chi-square
test_overall <- chisq.test(tab_intersp)
cat(sprintf("\nOverall chi-square: X2=%.2f, df=%d, P=%.6f\n",
            test_overall$statistic, test_overall$parameter, test_overall$p.value))

# Pairwise comparisons with Bonferroni correction (3 pairs)
# Use Fisher's exact for pairs involving L. dissortis (zero cell)
cat("\nPairwise comparisons (Bonferroni correction, n_pairs=3):\n")
species_names <- rownames(tab_intersp)
pairs <- list(c(1,2), c(1,3), c(2,3))

for (pr in pairs) {
  sp1 <- species_names[pr[1]]; sp2 <- species_names[pr[2]]
  sub  <- tab_intersp[pr, ]

  # Use Fisher's exact if any cell is 0, otherwise chi-square
  if (any(sub == 0)) {
    ft    <- fisher.test(sub)
    pval  <- ft$p.value
    stat_str <- sprintf("Fisher OR=%.2f", ft$estimate)
  } else {
    ct    <- chisq.test(sub)
    pval  <- ct$p.value
    stat_str <- sprintf("X2=%.2f", ct$statistic)
  }

  p_adj <- min(pval * 3, 1.0)
  sig   <- ifelse(p_adj < 0.001, "***",
            ifelse(p_adj < 0.01,  "**",
            ifelse(p_adj < 0.05,  "*", "ns")))
  cat(sprintf("  %s vs %s: %s, P=%.6f (adj. P=%.4f) %s\n",
              sp1, sp2, stat_str, pval, p_adj, sig))
}
