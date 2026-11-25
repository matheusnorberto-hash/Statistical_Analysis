# =========================================================
#@title CELL 2: STATISTICAL ANALYSIS (LMM / GLMM) AND BASIC PLOT
# =========================================================


library(dplyr)
library(ggplot2)
library(emmeans)
library(performance)
library(lme4)
library(lmerTest)

# 1. Filter ONLY Normal (excluding Log-Normal) and Gamma from the 'modelos' list
opcoes_validas <- names(modelos)[grepl("Normal|Gamma", names(modelos), ignore.case = TRUE) &
                                 !grepl("Log-Normal", names(modelos), ignore.case = TRUE)]

if (length(opcoes_validas) == 0) {
  stop("No valid models found. Please run Cell 1 first.")
}

cat("Available models for analysis:\n")
for (i in seq_along(opcoes_validas)) {
  cat(sprintf("[%d] %s\n", i, opcoes_validas[i]))
}

esc_dist <- as.integer(readline(prompt = "Enter the NUMBER of the selected distribution: "))
if (is.na(esc_dist) || esc_dist < 1 || esc_dist > length(opcoes_validas)) {
  stop(sprintf("Invalid choice! Enter a number between 1 and %d.", length(opcoes_validas)))
}
dist_escolhida <- opcoes_validas[esc_dist]

# Automatic identification of model family and link function
is_normal   <- grepl("Normal", dist_escolhida, ignore.case = TRUE) && !grepl("Log", dist_escolhida, ignore.case = TRUE)
is_log      <- grepl("log", dist_escolhida, ignore.case = TRUE)
tipo_modelo <- ifelse(is_normal, "LMM (Normal - Identity Link)", paste0("GLMM (", dist_escolhida, ")"))

# Check whether Cell 1 applied a Gamma transformation and flag it for this cell and downstream plots
transf_ativa <- (!is_normal) &&
                exists("transformacao_gamma_info") &&
                isTRUE(transformacao_gamma_info$aplicada)

cat(sprintf("\n>>> Selected Analysis: %s <<<\n", tipo_modelo))
if (transf_ativa) {
  cat("-------------------------------------------------------------------------\n")
  cat(sprintf("ℹ️ TRANSFORMATION FLAG DETECTED FROM CELL 1 (%s):\n", toupper(transformacao_gamma_info$tipo)))
  cat(sprintf("   %s\n", transformacao_gamma_info$descricao))
  cat("   Results below are reported on the original scale, and the flag remains active for figures.\n")
  cat("-------------------------------------------------------------------------\n")
}
cat("\n")

# 2. Final model fitting and Omnibus Test
if (is_normal) {
  modelo_final <- lmer(Valor ~ Condição * Momento + (1 | ID), data = dados_analise, REML = TRUE)
  omnibus <- anova(modelo_final, type = 3)

  # Total standard deviation for Cohen's d in LMM (Westfall et al., 2014)
  vc <- as.data.frame(VarCorr(modelo_final))
  sigma_total <- sqrt(sum(vc$vcov))
  df_resid <- df.residual(modelo_final)
} else {
  modelo_final <- modelos[[dist_escolhida]]
  omnibus <- emmeans::joint_tests(modelo_final)
}

# 3. Global Model Effect Size (Nakagawa's Marginal and Conditional R²)
r2_modelo <- tryCatch({
  res <- performance::r2(modelo_final)
  if (is.null(res) || is.na(res$R2_marginal)) NULL else res
}, error = function(e) NULL)

# Helper function to back-transform a numeric vector to the original scale when Gamma was transformed
reverter_escala_gamma <- function(x) {
  if (!transf_ativa) return(x)
  if (transformacao_gamma_info$tipo == "inversao") {
    return(x * -1)
  } else if (transformacao_gamma_info$tipo == "deslocamento") {
    return(x - transformacao_gamma_info$offset)
  }
  return(x)
}

# 4. Post-Hoc function with standardized Effect Size formatted to 2 decimal places
calcular_posthoc_completo <- function(emm_obj, ajuste = "bonferroni") {

  if (is_normal) {
    pares <- pairs(emm_obj, adjust = ajuste)
    df_res <- as.data.frame(pares)

    es <- emmeans::eff_size(emm_obj, sigma = sigma_total, edf = df_resid, method = "pairwise") %>%
      as.data.frame()
    lcl_col <- intersect(c("lower.CL", "asymp.LCL"), names(es))[1]
    ucl_col <- intersect(c("upper.CL", "asymp.UCL"), names(es))[1]

    df_res$Tamanho_Efeito <- sprintf("d = %.2f", es$effect.size)
    df_res$IC95_ES        <- sprintf("[%.2f, %.2f]", es[[lcl_col]], es[[ucl_col]])
    df_res$Delta_Pct      <- "-"

  } else if (is_log) {
    pares <- pairs(emm_obj, adjust = ajuste, type = "response")
    df_res <- as.data.frame(pares)
    names(df_res)[names(df_res) == "ratio"] <- "estimate"

    ci <- confint(pares, adjust = ajuste) %>% as.data.frame()
    lcl_col <- intersect(c("lower.CL", "asymp.LCL"), names(ci))[1]
    ucl_col <- intersect(c("upper.CL", "asymp.UCL"), names(ci))[1]

    df_res$Tamanho_Efeito <- sprintf("MR = %.2f", df_res$estimate)
    df_res$IC95_ES        <- sprintf("[%.2f, %.2f]", ci[[lcl_col]], ci[[ucl_col]])

    # If an additive shift (+constant) was applied, relative % change on the shifted scale is not applicable
    if (transf_ativa && transformacao_gamma_info$tipo == "deslocamento") {
      df_res$Delta_Pct <- "-"
    } else {
      df_res$Delta_Pct <- sprintf("%+.1f%%", (df_res$estimate - 1) * 100)
    }

  } else {
    # Gamma with Identity Link
    pares <- pairs(emm_obj, adjust = ajuste, type = "response")
    df_res <- as.data.frame(pares)

    ci <- confint(pares, adjust = ajuste) %>% as.data.frame()
    lcl_col <- intersect(c("lower.CL", "asymp.LCL"), names(ci))[1]
    ucl_col <- intersect(c("upper.CL", "asymp.UCL"), names(ci))[1]

    # If sign inversion (* -1) was applied, back-transform the difference and 95% CI to the original scale
    if (transf_ativa && transformacao_gamma_info$tipo == "inversao") {
      df_res$estimate <- df_res$estimate * -1
      lcl_orig <- ci[[ucl_col]] * -1
      ucl_orig <- ci[[lcl_col]] * -1
    } else {
      lcl_orig <- ci[[lcl_col]]
      ucl_orig <- ci[[ucl_col]]
    }

    df_res$Tamanho_Efeito <- sprintf("Diff = %.2f", df_res$estimate)
    df_res$IC95_ES        <- sprintf("[%.2f, %.2f]", lcl_orig, ucl_orig)

    emm_df <- as.data.frame(emm_obj, type = "response")
    col_val <- intersect(c("response", "emmean"), names(emm_df))[1]
    emm_df$val_orig <- reverter_escala_gamma(emm_df[[col_val]])

    grupos_comp <- strsplit(as.character(df_res$contrast), " - ")
    df_res$Delta_Pct <- sapply(seq_along(grupos_comp), function(idx) {
      if (transf_ativa && transformacao_gamma_info$tipo == "deslocamento") {
        return("-")
      }
      g2 <- trimws(grupos_comp[[idx]][2])
      if ("Momento" %in% names(df_res)) {
        mom <- df_res$Momento[idx]
        ref_val <- emm_df$val_orig[emm_df$Condição == g2 & emm_df$Momento == mom][1]
      } else {
        ref_val <- emm_df$val_orig[emm_df$Condição == g2][1]
      }
      if (is.na(ref_val) || ref_val == 0) "-" else sprintf("%+.1f%%", (df_res$estimate[idx] / ref_val) * 100)
    })
  }

  return(df_res)
}

# 5. Standardized and color-coded console table printer
imprimir_tabela_colorida <- function(df_dados) {
  RED <- "\033[1;31m"
  RESET <- "\033[0m"

  cols_interesse <- c("Momento", "contrast", "estimate", "t.ratio", "z.ratio",
                      "p.value", "Tamanho_Efeito", "IC95_ES", "Delta_Pct")
  cols_presentes <- intersect(cols_interesse, names(df_dados))
  df_clean <- df_dados[, cols_presentes, drop = FALSE]

  nomes_amigaveis <- sapply(names(df_clean), function(n) {
    switch(n,
           Momento        = "Bout",
           contrast       = "Contrast",
           estimate       = ifelse(is_log && !is_normal, "Ratio (MR)", "Difference"),
           t.ratio        = "t",
           z.ratio        = "z",
           p.value        = "p (Bonf.)",
           Tamanho_Efeito = "Effect Size",
           IC95_ES        = "95% CI (ES)",
           Delta_Pct      = "Change (%)",
           n)
  })

  cabecalho <- paste(sprintf("%-15s", nomes_amigaveis), collapse = " | ")
  cat(cabecalho, "\n")
  cat(paste(rep("-", nchar(cabecalho)), collapse = ""), "\n")

  for (i in 1:nrow(df_clean)) {
    valores_linha <- sapply(names(df_clean), function(col_name) {
      x <- df_clean[i, col_name]
      if (col_name == "p.value") {
        p_num <- as.numeric(x)
        if (!is.na(p_num) && p_num < 0.001) "< 0.001" else sprintf("%.4f", p_num)
      } else if (is.numeric(x)) {
        sprintf("%.2f", x)
      } else {
        as.character(x)
      }
    })

    linha_str <- paste(sprintf("%-15s", valores_linha), collapse = " | ")
    p_val <- as.numeric(df_clean$p.value[i])
    if (!is.na(p_val) && p_val < 0.05) {
      cat(RED, linha_str, " <- SIG", RESET, "\n", sep = "")
    } else {
      cat(linha_str, "\n")
    }
  }
}

# =========================================================
# 6. PRINT STATISTICAL RESULTS
# =========================================================
cat("=========================================================\n")
cat("1. OMNIBUS TEST -", tipo_modelo, "\n")
cat("=========================================================\n")
print(omnibus)

if (!is.null(r2_modelo)) {
  cat("\nGlobal Model Effect Size (Nakagawa):\n")
  cat(sprintf("-> Marginal R² (Fixed Effects):   %.3f (%.1f%%)\n", r2_modelo$R2_marginal, r2_modelo$R2_marginal * 100))
  cat(sprintf("-> Conditional R² (Full Model):   %.3f (%.1f%%)\n", r2_modelo$R2_conditional, r2_modelo$R2_conditional * 100))
}

cat("\n=========================================================\n")
cat("2. POST-HOC: MAIN EFFECT OF CONDITION (OVERALL MEAN)\n")
cat("=========================================================\n")
emm_cond <- emmeans(modelo_final, ~ Condição)
imprimir_tabela_colorida(calcular_posthoc_completo(emm_cond))

cat("\n=========================================================\n")
cat("3. POST-HOC: INTERACTION (CONDITION WITHIN EACH BOUT)\n")
cat("=========================================================\n")
emm_interacao <- emmeans(modelo_final, ~ Condição | Momento)
imprimir_tabela_colorida(calcular_posthoc_completo(emm_interacao))

# =========================================================
# 7. SCIENTIFIC AXIS LABELS (NIRS + METABOLIC + CPET)
# =========================================================
obter_eixo_ingles <- function(var_nome) {
  clean_name <- gsub("^(Esforço_|Recuperação_)", "", var_nome)

  # Identify phase (Exercise or Recovery)
  stage_prefix <- case_when(
    grepl("^Esforço|\\.E$", var_nome) ~ "Exercise",
    grepl("^Recuperação|\\.R$", var_nome) ~ "Recovery",
    TRUE ~ "Exercise"
  )

  # --- AUTOMATIC PARSER FOR NIRS VARIABLES (TSI, O2Hb, HHb, THb) ---
  if (grepl("^(TSI|O2Hb|HHb|THb)_(Ex|Ctrl)_", clean_name, ignore.case = TRUE)) {

    # 1. Identify leg (Exercised vs. Control)
    leg_str <- ifelse(grepl("_Ex_", clean_name, ignore.case = TRUE), "Exercised~Leg", "Control~Leg")

    # 2. Identify statistical/kinetic metric
    metric_str <- case_when(
      grepl("_Média$", clean_name, ignore.case = TRUE)     ~ "Mean",
      grepl("_Máx$", clean_name, ignore.case = TRUE)       ~ "Maximum",
      grepl("_Mín$", clean_name, ignore.case = TRUE)       ~ "Minimum",
      grepl("_Amplitude$", clean_name, ignore.case = TRUE) ~ "Amplitude",
      grepl("_Integral$", clean_name, ignore.case = TRUE)  ~ "Integral",
      grepl("Slope_5s_Max", clean_name, ignore.case = TRUE) ~ "Peak~5-s~Slope",
      grepl("t12_desox", clean_name, ignore.case = TRUE)    ~ "Deoxygenation~Half-Time~(t[1/2])",
      grepl("t12_reox", clean_name, ignore.case = TRUE)     ~ "Reoxygenation~Half-Time~(t[1/2])",
      grepl("Monoexp_Tau|Tau$", clean_name, ignore.case = TRUE) ~ "Time~Constant~(tau)",
      grepl("Monoexp_TD|TD$", clean_name, ignore.case = TRUE)   ~ "Time~Delay~(TD)",
      grepl("t_Nadir_s", clean_name, ignore.case = TRUE)    ~ "Time~to~Nadir",
      TRUE ~ ""
    )

    # 3. Identify NIRS signal and its respective unit (R plotmath: mu*M for μM)
    is_integral <- grepl("_Integral$", clean_name, ignore.case = TRUE)
    is_time     <- grepl("t12_|Tau$|TD$|t_Nadir", clean_name, ignore.case = TRUE)
    is_slope    <- grepl("Slope", clean_name, ignore.case = TRUE)

    signal_unit <- case_when(
      grepl("^TSI_", clean_name, ignore.case = TRUE) & is_integral  ~ "TSI~Integral~('%'*'\u00B7'*s)",
      grepl("^TSI_", clean_name, ignore.case = TRUE) & is_time      ~ "TSI~(s)",
      grepl("^TSI_", clean_name, ignore.case = TRUE) & is_slope     ~ "TSI~('%'*'\u00B7'*s^-1)",
      grepl("^TSI_", clean_name, ignore.case = TRUE)                ~ "TSI~('%')",

      grepl("^O2Hb_", clean_name, ignore.case = TRUE) & is_integral ~ "O[2]*Hb~Integral~(mu*M*'\u00B7'*s)",
      grepl("^O2Hb_", clean_name, ignore.case = TRUE)               ~ "O[2]*Hb~(mu*M)",

      grepl("^HHb_", clean_name, ignore.case = TRUE) & is_integral  ~ "HHb~Integral~(mu*M*'\u00B7'*s)",
      grepl("^HHb_", clean_name, ignore.case = TRUE)                ~ "HHb~(mu*M)",

      grepl("^THb_", clean_name, ignore.case = TRUE) & is_integral  ~ "tHb~Integral~(mu*M*'\u00B7'*s)",
      grepl("^THb_", clean_name, ignore.case = TRUE)                ~ "tHb~(mu*M)",
      TRUE ~ paste0("'", clean_name, "'")
    )

    # Avoid duplicating the word "Integral" when already present in the unit
    if (is_integral) {
      y_label <- paste0(leg_str, "~-~", signal_unit)
    } else if (metric_str != "") {
      y_label <- paste0(leg_str, "~-~", metric_str, "~", signal_unit)
    } else {
      y_label <- paste0(leg_str, "~-~", signal_unit)
    }

    return(list(label = y_label, stage = stage_prefix))
  }

  # --- METABOLIC AND CPET VARIABLES ---
  y_label <- case_when(
    grepl("^Potência$", clean_name, ignore.case = TRUE) ~ "Mean~Power~Output~(W)",
    grepl("Contribuição aerobica\\.mL\\.Kg", clean_name, ignore.case = TRUE)   ~ "Relative~Aerobic~Contribution~(mL*'\u00B7'*kg^-1)",
    grepl("Contribuição anaerobica\\.mL\\.Kg", clean_name, ignore.case = TRUE) ~ "Relative~Anaerobic~Contribution~(mL~O[2]~eq*'\u00B7'*kg^-1)",
    grepl("^Part\\.Aer\\.", clean_name, ignore.case = TRUE) ~ "Aerobic~Contribution~('%')",
    grepl("^Part\\.Ana\\.", clean_name, ignore.case = TRUE) ~ "Anaerobic~Contribution~('%')",
    grepl("^Contr\\.Aer\\.mL$", clean_name, ignore.case = TRUE) ~ "Absolute~Aerobic~Contribution~(mL)",
    grepl("^Contr\\.Ana\\.mL$", clean_name, ignore.case = TRUE) ~ "Absolute~Anaerobic~Contribution~(mL~O[2]~eq)",
    grepl("^Soma Aer\\+Ana$", clean_name, ignore.case = TRUE) ~ "Total~Energy~Expenditure~(mL~O[2]~eq)",
    grepl("^Sat\\.E$", clean_name, ignore.case = TRUE) ~ "Peripheral~Oxygen~Saturation~-~SpO[2]~('%')",
    grepl("^Sat\\.R$", clean_name, ignore.case = TRUE) ~ "Peripheral~Oxygen~Saturation~-~SpO[2]~('%')",
    grepl("^PSE$", clean_name, ignore.case = TRUE)   ~ "Rating~of~Perceived~Exertion~(CR-10)",

    grepl("VO2_Rel_Média", clean_name) ~ "Mean~VO[2]~(mL*'\u00B7'*kg^-1*'\u00B7'*min^-1)",
    grepl("VO2_Rel_Mín", clean_name)   ~ "Minimum~VO[2]~(mL*'\u00B7'*kg^-1*'\u00B7'*min^-1)",
    grepl("VO2_Rel_Máx", clean_name)   ~ "Maximum~VO[2]~(mL*'\u00B7'*kg^-1*'\u00B7'*min^-1)",
    grepl("VO2_Abs_Média", clean_name) ~ "Mean~Absolute~VO[2]~(L*'\u00B7'*min^-1)",
    grepl("VO2_Abs_Mín", clean_name)   ~ "Minimum~Absolute~VO[2]~(L*'\u00B7'*min^-1)",
    grepl("VO2_Abs_Máx", clean_name)   ~ "Maximum~Absolute~VO[2]~(L*'\u00B7'*min^-1)",
    grepl("VO2_Abs_Integral_L", clean_name)      ~ "Mean~VO[2]~Integral~(L)",
    grepl("VO2_Abs_Integral_mL_min", clean_name) ~ "VO[2]~Integral~(mL*'\u00B7'*min^-1)",
    grepl("FiO2_Média", clean_name) ~ "Mean~FiO[2]~('%')",
    grepl("FiO2_Mín", clean_name)   ~ "Minimum~FiO[2]~('%')",
    grepl("FiO2_Máx", clean_name)   ~ "Maximum~FiO[2]~('%')",
    grepl("VO2_HR", clean_name) ~ "Mean~Oxygen~Pulse~(mL*'\u00B7'*beat^-1)",
    grepl("VE_VO2_Média", clean_name) ~ "Mean~VE/VO[2]",
    grepl("VE_VO2_Mín", clean_name)   ~ "Minimum~VE/VO[2]",
    grepl("VE_VO2_Máx", clean_name)   ~ "Maximum~VE/VO[2]",
    grepl("VE_VCO2_Média", clean_name) ~ "Mean~VE/VCO[2]",
    grepl("VE_VCO2_Mín", clean_name)   ~ "Minimum~VE/VCO[2]",
    grepl("VE_VCO2_Máx", clean_name)   ~ "Maximum~VE/VCO[2]",
    grepl("HR_Média", clean_name) ~ "Mean~Heart~Rate~(bpm)",
    grepl("HR_Mín", clean_name)   ~ "Minimum~Heart~Rate~(bpm)",
    grepl("HR_Máx", clean_name)   ~ "Maximum~Heart~Rate~(bpm)",
    grepl("PetO2_Média", clean_name) ~ "Mean~PetO[2]~(mmHg)",
    grepl("PetO2_Mín", clean_name)   ~ "Minimum~PetO[2]~(mmHg)",
    grepl("PetO2_Máx", clean_name)   ~ "Maximum~PetO[2]~(mmHg)",
    grepl("PetCO2_Média", clean_name) ~ "Mean~PetCO[2]~(mmHg)",
    grepl("PetCO2_Mín", clean_name)   ~ "Minimum~PetCO[2]~(mmHg)",
    grepl("PetCO2_Máx", clean_name)   ~ "Maximum~PetCO[2]~(mmHg)",
    grepl("Duração_s", clean_name) ~ "Duration~(seconds)",
    TRUE ~ paste0("'", clean_name, "'")
  )
  return(list(label = y_label, stage = stage_prefix))
}

eixo_info <- obter_eixo_ingles(var_selecionada)

tabela_descritiva <- dados_analise %>%
  group_by(Condição, Momento) %>%
  summarise(
    n = n(),
    Média = mean(Valor, na.rm = TRUE),
    DP = sd(Valor, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Condição, as.numeric(as.character(Momento)))

p_basico <- ggplot(tabela_descritiva, aes(x = factor(Momento), y = Média, group = Condição, color = Condição)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Média - DP, ymax = Média + DP), width = 0.2) +
  labs(
    title = paste("Physiological Profile -", eixo_info$stage, "|", tipo_modelo),
    y = parse(text = eixo_info$label),
    x = "Bout (Effort)"
  ) +
  theme_light()

print(p_basico)
