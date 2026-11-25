# =========================================================
#@title CELL 1. LOAD DATA AND EVALUATE CANDIDATE DISTRIBUTIONS
# =========================================================


library(dplyr)
library(ggplot2)
library(glmmTMB)
library(patchwork)
library(knitr)

# 1. Interactive variable selection
escolha <- as.integer(readline(prompt = "Enter the NUMBER of the variable to analyze: "))
var_selecionada <- vars_validas[escolha]
cat("\n>>> Selected variable:", var_selecionada, "<<<\n\n")

# Prepare dataset removing NAs only for the selected variable
dados_analise <- df %>%
  select(ID, Condição, Momento, Valor = all_of(var_selecionada)) %>%
  filter(!is.na(Valor))

# ==============================================================================
# AUTOMATIC DETECTION AND TRANSFORMATION FOR GAMMA DISTRIBUTION
# ==============================================================================
val_min <- min(dados_analise$Valor, na.rm = TRUE)
val_max <- max(dados_analise$Valor, na.rm = TRUE)

if (val_max < 0) {
  # CASE 1: All values are strictly negative -> Sign inversion (* -1)
  dados_analise$Valor_Gamma <- dados_analise$Valor * -1
  transformacao_gamma_info <- list(
    aplicada  = TRUE,
    tipo      = "inversao",
    sinal     = -1,
    offset    = 0,
    val_min   = val_min,
    val_max   = val_max,
    descricao = sprintf("Sign inversion (x -1) applied because all values were negative (Min = %.4f | Max = %.4f).", val_min, val_max)
  )
  cat("=========================================================================\n")
  cat("🔄 NOTICE (GAMMA - SIGN INVERSION):\n")
  cat(sprintf("   All values for this variable are negative (Min = %.4f | Max = %.4f).\n", val_min, val_max))
  cat("   All values were inverted (x -1) to allow Gamma model fitting.\n")
  cat("   (This transformation has been flagged for downstream analysis and plots).\n")
  cat("=========================================================================\n\n")

} else if (val_min <= 0) {
  # CASE 2: Mixed values (negative and zero/positive) -> Add twice the absolute minimum
  constante_soma <- ifelse(val_min < 0, 2 * abs(val_min), 1.0)
  dados_analise$Valor_Gamma <- dados_analise$Valor + constante_soma
  transformacao_gamma_info <- list(
    aplicada  = TRUE,
    tipo      = "deslocamento",
    sinal     = 1,
    offset    = constante_soma,
    val_min   = val_min,
    val_max   = val_max,
    descricao = sprintf("Additive shift (+%.4f = twice the negative extreme %.4f) applied due to mixed-sign values.", constante_soma, val_min)
  )
  cat("=========================================================================\n")
  cat("➕ NOTICE (GAMMA - POSITIVE SHIFT FOR MIXED-SIGN VALUES):\n")
  cat(sprintf("   Mixed-sign values detected (negative and zero/positive). Negative extreme = %.4f.\n", val_min))
  cat(sprintf("   Added constant +%.4f (twice the negative extreme) to allow Gamma model fitting.\n", constante_soma))
  cat("   (This transformation has been flagged for downstream analysis and plots).\n")
  cat("=========================================================================\n\n")

} else {
  # CASE 3: All values are already strictly positive (> 0)
  dados_analise$Valor_Gamma <- dados_analise$Valor
  transformacao_gamma_info <- list(
    aplicada  = FALSE,
    tipo      = "nenhuma",
    sinal     = 1,
    offset    = 0,
    val_min   = val_min,
    val_max   = val_max,
    descricao = "No transformation required (all values were already > 0)."
  )
}

# Helper function to fit glmmTMB with optimizer fallback and convergence check
ajustar_glmm_seguro <- function(formula_mod, data_mod, familia_mod) {
  mod <- tryCatch({
    glmmTMB(formula_mod, data = data_mod, family = familia_mod, REML = FALSE)
  }, error = function(e) NULL)

  if (is.null(mod) || !isTRUE(mod$sdr$pdHess)) {
    mod <- tryCatch({
      glmmTMB(
        formula_mod, data = data_mod, family = familia_mod, REML = FALSE,
        control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
      )
    }, error = function(e) NULL)
  }

  if (!is.null(mod) && isTRUE(mod$sdr$pdHess)) {
    return(mod)
  }
  return(NULL)
}

# 2. Fit candidate models (all fitted via ML for valid AIC/BIC comparison)
modelos <- list()

modelos[["Normal (LMM)"]] <- glmmTMB(
  Valor ~ Condição * Momento + (1 | ID),
  data = dados_analise,
  family = gaussian(link = "identity"),
  REML = FALSE
)

# A) Gamma with Log Link (fitted on Valor_Gamma)
mod_gamma_log <- ajustar_glmm_seguro(
  Valor_Gamma ~ Condição * Momento + (1 | ID),
  dados_analise,
  Gamma(link = "log")
)
if (!is.null(mod_gamma_log)) {
  modelos[["Gamma (Link Log)"]] <- mod_gamma_log
} else {
  cat("WARNING: Gamma (Link Log) did not converge stably.\n")
}

# B) Gamma with Identity Link (fitted on Valor_Gamma)
mod_gamma_id <- ajustar_glmm_seguro(
  Valor_Gamma ~ Condição * Momento + (1 | ID),
  dados_analise,
  Gamma(link = "identity")
)
if (!is.null(mod_gamma_id)) {
  modelos[["Gamma (Link Identity)"]] <- mod_gamma_id
} else {
  cat("WARNING: Gamma (Link Identity) did not converge stably.\n")
}

# 3. Comparative table of AIC, BIC, and Log-Likelihood
comparacao <- data.frame(
  Distribuição = character(),
  AIC = numeric(),
  BIC = numeric(),
  LogLik = numeric(),
  stringsAsFactors = FALSE
)

for (nome in names(modelos)) {
  m <- modelos[[nome]]
  comparacao <- rbind(comparacao, data.frame(
    Distribuição = nome,
    AIC = round(AIC(m), 2),
    BIC = round(BIC(m), 2),
    LogLik = round(as.numeric(logLik(m)), 2)
  ))
}

comparacao <- comparacao %>% arrange(AIC)
comparacao$Delta_AIC <- round(comparacao$AIC - min(comparacao$AIC), 2)

cat("=== MODEL FIT COMPARISON (Lower AIC and BIC = Better Fit) ===\n\n")
print(knitr::kable(
  comparacao,
  format = "simple",
  align = "llrrrr",
  col.names = c("Model / Distribution", "AIC", "BIC", "LogLik", "Delta AIC")
))
cat("\nRule of thumb: Delta AIC > 2 indicates substantial support for the model with lower AIC.\n\n")

# 4. Residual Diagnostic Plots (Fitted values back-transformed to original scale)
lista_plots <- list()

for (nome in names(modelos)) {
  m <- modelos[[nome]]
  res <- residuals(m, type = "pearson")
  pred_bruto <- fitted(m)

  # Back-transform fitted values on the X-axis if a Gamma transformation was applied
  if (grepl("Gamma", nome) && isTRUE(transformacao_gamma_info$aplicada)) {
    if (transformacao_gamma_info$tipo == "inversao") {
      pred <- pred_bruto * -1
    } else if (transformacao_gamma_info$tipo == "deslocamento") {
      pred <- pred_bruto - transformacao_gamma_info$offset
    } else {
      pred <- pred_bruto
    }
  } else {
    pred <- pred_bruto
  }

  df_diag <- data.frame(Residuos = res, Preditos = pred)

  p1 <- ggplot(df_diag, aes(x = Residuos)) +
    geom_histogram(aes(y = after_stat(density)), bins = 25, fill = "#4E79A7", alpha = 0.6, color = "white") +
    geom_density(color = "#E15759", linewidth = 1.2) +
    stat_function(fun = dnorm, args = list(mean = mean(res), sd = sd(res)), linetype = "dashed", color = "black") +
    labs(title = paste(nome, "- Residual Density"), x = "Pearson Residuals", y = "Density") +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(size = 13, face = "bold"), plot.margin = margin(10, 15, 10, 15))

  p2 <- ggplot(df_diag, aes(sample = Residuos)) +
    stat_qq(color = "#4E79A7", alpha = 0.7, size = 1.8) +
    stat_qq_line(color = "#E15759", linewidth = 1.2) +
    labs(title = paste(nome, "- Q-Q Plot"), x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(size = 13, face = "bold"), plot.margin = margin(10, 15, 10, 15))

  p3 <- ggplot(df_diag, aes(x = Preditos, y = Residuos)) +
    geom_point(color = "#4E79A7", alpha = 0.6, size = 1.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "#E15759", linewidth = 1.2) +
    geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 1) +
    labs(title = paste(nome, "- Fitted vs. Residuals"), x = "Fitted Values (Original Scale)", y = "Pearson Residuals") +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(size = 13, face = "bold"), plot.margin = margin(10, 15, 10, 15))

  p_unido <- (p1 | p2 | p3) + plot_layout(ncol = 3, widths = c(1, 1, 1))
  lista_plots[[length(lista_plots) + 1]] <- p_unido
}

print(wrap_plots(lista_plots, ncol = 1))
