################################################################################
# Correation analyses
################################################################################
#
# Lina Aragón
# linamaragonb@gmail.com
# September 24, 2025
#
# Description
################################################################################

# 1. Load libraries ----

#install.packages("devtools", "tidyverse", "ggplot2", "plotly", "vegan", "ggvegan", "ggrepel", "glue", "viridis", "fixest", "lmtest", "corrplot") #install if needed
# devtools::install_github("gavinsimpson/ggvegan")
pkgs <- c("devtools", "tidyverse", "ggplot2", "plotly", "vegan", "ggvegan", "ggrepel", "glue", "viridis", "fixest",
          "lmtest", "corrplot")
lapply(pkgs, library, character.only = TRUE)
remove(pkgs)

# 2. Load data ----
trait_data <- read_csv("data/processed/v_PFCT7_clean_functional_traits_2023.csv")
trait_data_wide <- read_csv("data/processed/v_PFCT7_clean_functional_traits_2023_wide.csv")



# 3. Trait Correlation ----

## Matrix of traits ----
matrix_traits <- trait_data_wide |> 
  select(root_depth, veg_height, rd, bi, srl, rtd, rdmc, sla, ldmc, bgb_agb) |>
  mutate(across(everything(), as.numeric))

## Correlation matrix ----
cor_matrix_traits <- cor(matrix_traits, use = "pairwise.complete.obs", method = "pearson")

## Compute p-value matrix ----
cor_p_matrix <- function(x) {
  n <- ncol(x)
  pmat <- matrix(NA, n, n)
  colnames(pmat) <- rownames(pmat) <- colnames (x)
  for (i in 1:n) {
    for(j in 1:n) {
      if (i == j) {
        pmat[i,j] <- NA
      } else {
        pmat[i,j] <- cor.test(x[[i]],x[[j]])$p.value
      }
    }
  }
  pmat
}

pval_matrix_traits <- cor_p_matrix(matrix_traits)

# Diagnostics
#stopifnot(
#  identical(dim(cor_matrix_traits), dim(pval_matrix_traits)),
#  all(rownames(cor_matrix_traits) == rownames(pval_matrix_traits)),
#  all(colnames(cor_matrix_traits) == colnames(pval_matrix_traits))
#)
#
#diag(cor_matrix_traits) <- 1
#diag(pval_matrix_traits) <- NA

# Visual representation
corrplot(
  cor_matrix_traits,
  method = "color",
  type = "lower",
  tl.srt = 45,
  addCoef.col = "black",      # Show coefficients
  tl.col = "black",           # Axis text color
  p.mat = pval_matrix_traits, # P-value matrix
  sig.level = 0.05,           # Significance threshold
  insig = "blank",            # Only blanks numbers, keeps colors for all
  number.digits = 2           # Number of decimals for coefficients
)

# 4. Root trait relationships across species using linear, polynomial and exponential model ----

# Trait pairs
pairs <- tribble(
  ~trait1, ~trait2,
  "rd", "srl",
  "root_depth", "rtd",
  "bi", "rtd",
  "rtd", 'srl', 
  "rdmc", "rtd",
  "root_depth", "sla",
  "bi", "sla",
  "root_depth", "ldmc",
  "bi", "ldmc",
  "ldmc", "sla",
  "bgb_agb", "ldmc"
)


# Model types
formulas <- tribble(
  ~type, ~formula_template,
  "linear", "{trait2} ~ {trait1}",
  "poly", "{trait2} ~ {trait1} + I({trait1}^2)",
  "exp", "log({trait2}) ~ {trait1}"
)

# Expand for all pairs and formulas
model_grid <- crossing(pairs, formulas) %>%
  mutate(
    base_formula = pmap_chr(
      list(formula_template, trait1, trait2),
      ~ glue(.x, trait1 = .y, trait2 = ..3)
    ),
    model_name = paste(trait1, trait2, type, sep = "_")
  )

# Fit all models
results <- model_grid %>%
  mutate(
    model = map(base_formula, ~feols(as.formula(.x), data = trait_data_wide)),
    AIC = map_dbl(model, AIC),
    BIC = map_dbl(model, BIC),
    summary = map(model, summary),
    residuals = map(model, resid),
    fitted = map(model, fitted),
    shapiro = map(residuals, ~shapiro.test(.x)),
    shapiro_W = map_dbl(shapiro, ~ if(is.null(.x) || is.null(.x$statistic)) NA_real_ else as.numeric(.x$statistic)),
    bp_test = map2(residuals, fitted, ~bptest(.x ~ .y)),
    r_squared = map_dbl(model, ~fixest::r2(.x, type = "r2"))
  )


# Summarize results

mod <- feols(srl ~ rtd, trait_data_wide)
summary(mod)
summary(mod, cluster =  ~species)

mod$coeftable


summary_table <- model_grid %>%
  mutate(
    model = map(base_formula, ~feols(as.formula(.x), data = trait_data_wide)),
    summary_obj = map(model, ~ if (inherits(.x, "try-error")) NA else summary(.x)),
    # Extract non-intercept term estimates, p-values, and significance
    coef_df = map(summary_obj, function(s) {
      if (is.null(s) || all(is.na(s))) {
        tibble(term = NA_character_, estimate = NA_real_, p.value = NA_real_, significant = NA)
      } else {
        ct <- as.data.frame(s$coeftable)
        ct$term <- rownames(ct)
        ct %>%
          filter(term != "(Intercept)") %>%
          select(term, estimate = Estimate, p.value = `Pr(>|t|)`) %>%
          mutate(significant = if_else(!is.na(p.value) & p.value < 0.05, TRUE, FALSE))
      }
    }),
    AIC = map_dbl(model, AIC),
    BIC = map_dbl(model, BIC),
    residuals = map(model, resid),
    fitted = map(model, fitted),
    shapiro = map(residuals, ~shapiro.test(.x)),
    bp_test = map2(residuals, fitted, ~bptest(.x ~ .y)),
    r_squared = map_dbl(model, ~if (inherits(.x, "try-error")) NA_real_ else fixest::r2(.x, type = "r2"))
  ) %>%
  transmute(
    model_name,
    formula = base_formula,
    warning = map_chr(model, ~ifelse(any(class(.x) == "try-error"), "ERROR", "")),
    coef_df,
    AIC,
    BIC,
    shapiro_W = map_dbl(shapiro, ~ if(is.null(.x) || is.null(.x$statistic)) NA_real_ else as.numeric(.x$statistic)),
    shapiro_p = map_dbl(shapiro, "p.value"),
    shapiro_pass = shapiro_p > 0.05,
    bp_p = map_dbl(bp_test, "p.value"),
    bp_pass = bp_p > 0.05,
    r_squared
  ) %>%
  unnest(coef_df)

print(summary_table)

summary_table <- model_grid %>%
  mutate(
    # Add species as a fixed effect (dummy variable coding)
    model = map(base_formula, ~feols(as.formula(.x),, data = trait_data_wide)),
    summary_obj = map(model, ~if (inherits(.x, "try-error")) NA else summary(.x)),
    # Calculate overall F-statistic and p-value
    F_value = map_dbl(summary_obj, function(s) {
      if (is.null(s) || all(is.na(s))) return(NA_real_)
      nparams <- s$nparams
      nobs <- s$nobs
      ssr_null <- s$ssr_null
      ssr <- s$ssr
      df1 <- nparams - 1
      df2 <- nobs - nparams
      num <- (ssr_null - ssr) / df1
      den <- ssr / df2
      Fstat <- num / den
      if (is.na(Fstat) || is.nan(Fstat) || is.infinite(Fstat)) return(NA_real_)
      Fstat
    }),
    F_p = map2_dbl(summary_obj, F_value, function(s, Fval) {
      if (is.null(s) || all(is.na(s)) || is.na(Fval)) return(NA_real_)
      nparams <- s$nparams
      nobs <- s$nobs
      df1 <- nparams - 1
      df2 <- nobs - nparams
      pf(Fval, df1, df2, lower.tail = FALSE)
    }),
    model_significant = F_p < 0.05,
    r_squared = map_dbl(model, ~if (inherits(.x, "try-error")) NA_real_ else fixest::r2(.x, type = "r2")),
    AIC = map_dbl(model, ~if (inherits(.x, "try-error")) NA_real_ else AIC(.x)),
    BIC = map_dbl(model, ~if (inherits(.x, "try-error")) NA_real_ else BIC(.x))
  ) %>%
  transmute(
    model_name,
    formula = base_formula,
    warning = map_chr(model, ~ifelse(any(class(.x) == "try-error"), "ERROR", "")),
    F_value,
    F_p,
    model_significant,
    r_squared,
    AIC,
    BIC
  )

summary_table_spp <- model_grid %>%
  mutate(
    # Add species as a fixed effect (dummy variable coding)
    model = map(base_formula, ~feols(update(as.formula(.x), . ~ . + species), data = trait_data_wide)),
    summary_obj = map(model, ~if (inherits(.x, "try-error")) NA else summary(.x)),
    # Calculate overall F-statistic and p-value
    F_value = map_dbl(summary_obj, function(s) {
      if (is.null(s) || all(is.na(s))) return(NA_real_)
      nparams <- s$nparams
      nobs <- s$nobs
      ssr_null <- s$ssr_null
      ssr <- s$ssr
      df1 <- nparams - 1
      df2 <- nobs - nparams
      num <- (ssr_null - ssr) / df1
      den <- ssr / df2
      Fstat <- num / den
      if (is.na(Fstat) || is.nan(Fstat) || is.infinite(Fstat)) return(NA_real_)
      Fstat
    }),
    F_p = map2_dbl(summary_obj, F_value, function(s, Fval) {
      if (is.null(s) || all(is.na(s)) || is.na(Fval)) return(NA_real_)
      nparams <- s$nparams
      nobs <- s$nobs
      df1 <- nparams - 1
      df2 <- nobs - nparams
      pf(Fval, df1, df2, lower.tail = FALSE)
    }),
    model_significant = F_p < 0.05,
    r_squared = map_dbl(model, ~if (inherits(.x, "try-error")) NA_real_ else fixest::r2(.x, type = "r2")),
    AIC = map_dbl(model, ~if (inherits(.x, "try-error")) NA_real_ else AIC(.x)),
    BIC = map_dbl(model, ~if (inherits(.x, "try-error")) NA_real_ else BIC(.x))
  ) %>%
  transmute(
    model_name,
    formula = base_formula,
    warning = map_chr(model, ~ifelse(any(class(.x) == "try-error"), "ERROR", "")),
    F_value,
    F_p,
    model_significant,
    r_squared,
    AIC,
    BIC
  )

## Model variations
#variations <- tribble(
#  ~variation, ~suffix,
#  "none", "",
#  "species_FE", "| species",
#  "species_interact", "* species"
#)

## Expand for all pairs, formulas, and variations
#model_grid <- crossing(pairs, formulas, variations) %>%
#  mutate(
#    base_formula = pmap_chr(
#      list(formula_template, trait1, trait2),
#      ~ glue(.x, trait1 = .y, trait2 = ..3)
#    ),
#    full_formula = case_when(
#      variation == "species_interact" ~ 
#        paste0(sub(" ~ ", " ~ (", base_formula, fixed = TRUE), ") * species"),
#      TRUE ~ paste(base_formula, suffix)
#    ),
#    model_name = paste(trait1, trait2, type, variation, sep = "_")
#  )

## Fit all models
#results <- model_grid %>%
# mutate(
#   model = map(full_formula, ~feols(as.formula(.x), data = trait_data_wide)),
#   AIC = map_dbl(model, AIC),
#   BIC = map_dbl(model, BIC),
#   summary = map(model, summary),
#   residuals = map(model, resid),
#   fitted = map(model, fitted),
#   shapiro = map(residuals, ~shapiro.test(.x)),
#   bp_test = map2(residuals, fitted, ~bptest(.x ~ .y)),
#   r_squared = map_dbl(model, ~fixest::r2(.x, type = "r2"))
# )


# 4. Trait Correlation ----

#Matrix of traits
matrix_traits <- trait_data_wide |> 
  select(root_depth, veg_height, rd, bi, srl, rtd, rdmc, sla, ldmc, bgb_agb) |>
  mutate(across(everything(), as.numeric))

#Correlation matrix
cor_matrix_traits <- cor(matrix_traits, use = "pairwise.complete.obs", method = "pearson")

#Compute p-value matrix
cor_p_matrix <- function(x) {
  n <- ncol(x)
  pmat <- matrix(NA, n, n)
  colnames(pmat) <- rownames(pmat) <- colnames (x)
  for (i in 1:n) {
    for(j in 1:n) {
      if (i == j) {
        pmat[i,j] <- NA
      } else {
        pmat[i,j] <- cor.test(x[[i]],x[[j]])$p.value
      }
    }
  }
  pmat
}

pval_matrix_traits <- cor_p_matrix(matrix_traits)

# Diagnostics
stopifnot(
  identical(dim(cor_matrix_traits), dim(pval_matrix_traits)),
  all(rownames(cor_matrix_traits) == rownames(pval_matrix_traits)),
  all(colnames(cor_matrix_traits) == colnames(pval_matrix_traits))
)

diag(cor_matrix_traits) <- 1
diag(pval_matrix_traits) <- NA


# Visual representation
corrplot(
  cor_matrix_traits,
  method = "color",
  type = "lower",
  tl.srt = 45,
  addCoef.col = "black",      # Show coefficients
  tl.col = "black",           # Axis text color
  p.mat = pval_matrix_traits, # P-value matrix
  sig.level = 0.05,           # Significance threshold
  insig = "blank",            # Only blanks numbers, keeps colors for all
  number.digits = 2           # Number of decimals for coefficients
)


# 4. Similar analysis to Weemstra et al 2020 Functional ecology ----

# Your trait pairs
pairs <- tribble(
  ~trait1, ~trait2,
  "rd", "srl",
  "root_depth", "rtd",
  "bi", "rtd",
  "rtd", 'srl', 
  "rdmc", "rtd",
  "root_depth", "sla",
  "bi", "sla",
  "root_depth", "ldmc",
  "bi", "ldmc",
  "ldmc", "sla",
  "bgb_agb", "ldmc"
)

# Model types
formulas <- tribble(
  ~type, ~formula_template,
  "linear", "{trait2} ~ {trait1}",
  "poly", "{trait2} ~ {trait1} + I({trait1}^2)",
  "exp", "log({trait2}) ~ {trait1}"
)

# Model variations
variations <- tribble(
  ~variation, ~suffix,
  "none", "",
  "species_FE", "| species",
  "species_interact", "* species"
)

# Expand for all pairs, formulas, and variations
model_grid <- crossing(pairs, formulas, variations) %>%
  mutate(
    base_formula = pmap_chr(
      list(formula_template, trait1, trait2),
      ~ glue(.x, trait1 = .y, trait2 = ..3)
    ),
    full_formula = case_when(
      variation == "species_interact" ~ 
        paste0(sub(" ~ ", " ~ (", base_formula, fixed = TRUE), ") * species"),
      TRUE ~ paste(base_formula, suffix)
    ),
    model_name = paste(trait1, trait2, type, variation, sep = "_")
  )

# Fit all models
results <- model_grid %>%
  mutate(
    model = map(full_formula, ~feols(as.formula(.x), data = trait_data_wide)),
    AIC = map_dbl(model, AIC),
    BIC = map_dbl(model, BIC),
    summary = map(model, summary),
    residuals = map(model, resid),
    fitted = map(model, fitted),
    shapiro = map(residuals, ~shapiro.test(.x)),
    bp_test = map2(residuals, fitted, ~bptest(.x ~ .y)),
    r_squared = map_dbl(model, ~fixest::r2(.x, type = "r2"))
  )

# Summarize results
summary_table <- results %>%
  transmute(
    model_name,
    formula = full_formula,
    warning = map_chr(model, ~ifelse(any(class(.x) == "try-error"), "ERROR", "")),
    shapiro_p = map_dbl(shapiro, "p.value"),
    shapiro_pass = shapiro_p > 0.05,
    bp_p = map_dbl(bp_test, "p.value"),
    bp_pass = bp_p > 0.05,
    r_squared
  )
print(summary_table)

# 1. Collect all AIC/BIC values for all models
aic_bic_table <- results %>%
  select(trait1, trait2, model_name, full_formula, AIC, BIC)

# 2. For each trait pair, get the model with the minimum AIC and BIC
best_models <- aic_bic_table %>%
  group_by(trait1, trait2) %>%
  summarise(
    best_AIC = min(AIC, na.rm = TRUE),
    best_AIC_model = model_name[which.min(AIC)],
    best_BIC = min(BIC, na.rm = TRUE),
    best_BIC_model = model_name[which.min(BIC)],
    .groups = "drop"
  )

# 3. If you want a single wide table combining all info:
final_table <- aic_bic_table %>%
  left_join(best_models, by = c("trait1", "trait2"))

# 4. If you want just the best model rows per trait pair, do:
best_only <- best_models

# Print the best model for each trait pair (by AIC and BIC)
print(best_only)



results %>%
  group_by(trait1, trait2) %>%
  mutate(
    aic_rank = rank(AIC, ties.method = "min"),
    bic_rank = rank(BIC, ties.method = "min"),
    combined_rank = aic_rank + bic_rank
  ) %>%
  filter(combined_rank == min(combined_rank, na.rm = TRUE)) %>%
  ungroup() %>%
  select(trait1, trait2, model_name, full_formula, AIC, BIC, combined_rank)
