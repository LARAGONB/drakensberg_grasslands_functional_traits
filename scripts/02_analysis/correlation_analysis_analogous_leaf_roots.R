################################################################################
# Leaf vs Roots Trait correlation across all plants and within species
################################################################################
#
# Lina Aragón
# linamaragonb@gmail.com
# September 24, 2025
#
# Description
# Here we evaluated trait correlations across all the plants measured for the
# 5 species of interest and within it species. We used AIC to find the best model
# linear, polynomial or exponential, and create a graph including all the traits
# combinations that showed to be significant the wether or not they indeed had a 
# significant relationship across all plants and within species
################################################################################

# 1. Load libraries ----

# install.packages("devtools", "tidyverse", "ggplot2", "plotly", "vegan", "ggvegan", "ggrepel",
#                  "glue", "viridis", "fixest", "lmtest", "corrplot", "skedastic","patchwork",
#                  "cowplot") #install if needed
# devtools::install_github("gavinsimpson/ggvegan")
pkgs <- c("devtools", "tidyverse", "ggplot2", "plotly", "vegan", "ggvegan", "ggrepel",
          "glue", "viridis", "fixest", "lmtest", "corrplot", "skedastic", "patchwork",
          "cowplot")
lapply(pkgs, library, character.only = TRUE)
remove(pkgs)

# 2. Load data ----
trait_data <- read_csv("data/processed/v_PFCT7_clean_functional_traits_2023.csv")
trait_data_wide <- read_csv("data/processed/v_PFCT7_clean_functional_traits_2023_wide.csv")

# 3. Define traits ----
traits <- c("leaf_thickness", "sla", "ldmc", "rd", "srl", "rtd")

metadata <- c("id", "aspect", "site_id", "elevation_m_asl", "plant_id", "species", "family", "growth_form")

# 4. Generate ALL possible pairs----
pairs <- expand_grid(trait_x = traits, trait_y = traits) |> 
  filter(trait_x != trait_y) |>  # Remove self-pairs
  mutate(traits = paste(trait_y, trait_x, sep = "~")) |> 
  filter(traits %in% c("sla~srl", "srl~sla", 
                       "rd~leaf_thickness", "leaf_thickness~rd",
                       "rtd~ldmc", "ldmc~rtd"))

# 5. Model types ----
formulas <- tribble(
  ~type, ~formula_template,
  "linear", "{trait_y} ~ {trait_x}",
  "poly", "{trait_y} ~ {trait_x} + I({trait_x}^2)",
  "exp", "log({trait_y} + 1e-10) ~ {trait_x}"
)

# 5. Expand for all pairs and formulas ----
model_grid <- crossing(pairs, formulas) |> 
  mutate(
    base_formula = pmap_chr(
      list(formula_template, trait_x, trait_y),
      ~ glue(.x, trait_x = .y, trait_y = ..3)
    ),
    model_name = paste(traits, type, sep = "_")
  ) |> 
  select(!formula_template)

model_grid

# 6. Fit all models and extract valuable information ----
results <- model_grid |> 
  mutate(
    model = map(base_formula, ~ feols(as.formula(.x), data = trait_data_wide)),
    estimate = map_dbl(model, ~ coef(.x)[[2]], se = "hc3",),
    mean_dv = map_dbl(model, ~ fitstat(.x, "my")$my, se = "hc3",),
    # se = map_dbl(model, ~ fixest::se(.x)[[2]]),
    rmse = map_dbl(model, ~ fitstat(.x, "rmse")$rmse, se = "hc3",),
    f_value = map_dbl(model, ~ fitstat(.x, "f")$f$stat, se = "hc3",),
    p_value = map_dbl(model, ~ fitstat(.x, "f")$f$p, se = "hc3",),
    ar2 = map_dbl(model, ~ fitstat(.x, "ar2")$ar2, se = "hc3",),
    significant = p_value < 0.05,
    aic = map_dbl(model, ~ fitstat(.x, "aic")$aic, se = "hc3",),
    w_shapiro = map_dbl(model, ~ shapiro.test(resid(.x))$statistic[[1]]),
    p_shapiro = map_dbl(model, ~ shapiro.test(resid(.x))$p.value),
    w_shap_normal = w_shapiro > 0.9,
    p_shap_normal = p_shapiro >= 0.05,
    s_white = map_dbl(base_formula, ~ white(lm(as.formula(.x), data = trait_data_wide))$statistic),
    p_white = map_dbl(base_formula, ~ white(lm(as.formula(.x), data = trait_data_wide))$p.value),
    p_whit_homosked = p_white >= 0.05,
    p_bp = map_dbl(model, ~ bptest(resid(.x) ~ fitted(.x))$p.value),
    p_bp_homosked = p_bp >= 0.05,
    etable = map(model, ~ etable(.x, se = "hc3", fitstat = ~ . + f + my + rmse)),
    predicted = map2(model, type, ~ 
                       if (.y %in% c("linear", "poly")) {
                         predict(.x)
                       } else if (.y == "exp") {
                         exp(predict(.x))
                       } else {
                         NA  # or handle other model types
                       })) |> 
  relocate(model, .before = etable)
results  

## Sort table based on aic value for each trait pair ----
results_sort <- results |> 
  group_by(traits) |> 
  arrange(traits,aic)

## Extract best model (lowest aic) for each trait pair ----
results_best <- results %>%
  group_by(traits) %>%
  filter(aic == min(aic)) %>%
  ungroup()

## Extract significant model ----

sig_best <- results_best |> 
  filter(significant == TRUE)

### Export table ----
sig_best |> 
  select(traits, type, estimate, f_value, p_value, ar2, w_shapiro, s_white, p_white, p_bp) |> 
  arrange(traits) |> 
  write_csv("results/tab/trait_relationship_across_leaf_roots.csv")

## Extract the fitted values per trait pair ----
selected_cols <- c(metadata, traits, "predicted")

results_blong <- results_best |> 
  mutate(
    data = map(predicted, ~ mutate(trait_data_wide, predicted = .x))
  )  |> 
  select(1:6, data) |> 
  unnest(data) |> 
  select(1:6, all_of(selected_cols)) |> 
  rename(predicted_all = predicted)

# 7. Within species correlations ----

## Nest data by species
nested_spp <- trait_data_wide  |> 
  group_by(species)  |> 
  nest()

## Created grid including species 
species_grid <- crossing(
  species = unique(trait_data_wide$species),
  model_grid) |> 
  left_join(nested_spp, by = "species")


## Fit all models and extract valuable information ----
results_spp <- species_grid |> 
  mutate(
    model = map2(base_formula, data, ~ feols(as.formula(.x), data = .y)),
    estimate = map_dbl(model, ~ coef(.x)[[2]], se = "hc3",),
    mean_dv = map_dbl(model, ~ fitstat(.x, "my")$my, se = "hc3"),
    se = map_dbl(model, ~ fixest::se(.x)[[2]], se = "hc3"),
    rmse = map_dbl(model, ~ fitstat(.x, "rmse")$rmse, se = "hc3"),
    f_value = map_dbl(model, ~ fitstat(.x, "f")$f$stat, se = "hc3"),
    p_value = map_dbl(model, ~ fitstat(.x, "f")$f$p, se = "hc3"),
    ar2 = map_dbl(model, ~ fitstat(.x, "ar2")$ar2, se = "hc3"),
    significant = p_value < 0.05,
    aic = map_dbl(model, ~ fitstat(.x, "aic")$aic, se = "hc3"),
    w_shapiro = map_dbl(model, ~ shapiro.test(resid(.x))$statistic[[1]]),
    p_shapiro = map_dbl(model, ~ shapiro.test(resid(.x))$p.value),
    w_shap_normal = w_shapiro > 0.9,
    p_shap_normal = p_shapiro >= 0.05,
    s_white = map2_dbl(base_formula, data, ~ white(lm(as.formula(.x), data = .y))$statistic),
    p_white = map2_dbl(base_formula, data, ~ white(lm(as.formula(.x), data = .y))$p.value),
    p_whit_normal = p_white >= 0.05,
    p_bp = map_dbl(model, ~bptest(resid(.x) ~ fitted(.x))$p.value),
    p_bp_normal = p_bp >= 0.05,
    etable = map(model, ~ etable(.x, se = "hc3", fitstat = ~ . + f + my + rmse)),
    predicted = map2(model, type, ~ 
                       if (.y %in% c("linear", "poly")) {
                         predict(.x)
                       } else if (.y == "exp") {
                         exp(predict(.x))
                       } else {
                         NA  # or handle other model types
                       })) |> 
  relocate(model, .before = etable)
results_spp  

## Sort table based on aic value for each trait pair ----
results_sort_spp <- results_spp |> 
  group_by(species, traits) |> 
  arrange(species, traits, aic)

## Extract best model (lowest aic) for each trait pair ----
results_best_spp <- results_spp %>%
  group_by(species, traits) %>%
  filter(aic == min(aic)) %>%
  ungroup()

## Extract significant model ----

sig_best_spp <- results_best_spp |> 
  filter(significant == TRUE)

### Export table ----
sig_best_spp |> 
  select(species, traits, type, estimate, f_value, p_value, ar2, w_shapiro, s_white, p_white, p_bp) |> 
  arrange(species, traits) |> 
  write_csv("results/tab/trait_relationship_within_leaf_roots.csv")


## Extract the fitted values per species and trait pair ----
selected_cols_spp <- c(metadata, traits, "predicted")

results_blong_spp <- results_best_spp %>%
  mutate(
    data = map2(data, predicted, ~ mutate(.x, predicted = .y))
  ) %>%
  select(1:7, data) %>%   # First 7 cols + data
  unnest(data) |> 
  select(1:7, all_of(selected_cols_spp)) |> 
  rename(predicted_spp = predicted)

# 8. Full results table ----

join_cols <- setdiff(intersect(names(results_blong), names(results_blong_spp)), c("type", "base_formula", "model_name"))

results_general <- results_blong |> 
  left_join(
    results_blong_spp |> select(-any_of(c("type", "base_formula", "model_name"))),
    by = c("traits", "species", "id", "aspect", "site_id", "elevation_m_asl", 
           "plant_id", "family", "growth_form", traits)) |> 
  mutate(species = fct_relevel(species, "Themeda triandra", after = 2))


# 9. Trait pairs significant in EITHER analysis ----
sig_trait_pairs_general <- sig_best |> 
  pull(traits) |> 
  unique()

sig_trait_pairs_species <- sig_best_spp |> 
  pull(traits) |> 
  unique()

all_sig_trait_pairs <- c(sig_trait_pairs_general, sig_trait_pairs_species) |> 
  unique()

## Create classification of trait pairs ----
trait_pair_classification <- tibble(
  traits = all_sig_trait_pairs) |> 
  mutate(
    sig_general = traits %in% sig_trait_pairs_general,
    sig_species = traits %in% sig_trait_pairs_species,
    n_sig_species = map_int(traits, ~ {
      sig_best_spp |> filter(traits == .x) |> nrow()}),
    category = case_when(
      sig_general & sig_species ~ "Both",
      sig_general & !sig_species ~ "General only",
      !sig_general & sig_species ~ "Species only",
      TRUE ~ "Neither"))

## View summary
trait_pair_classification |> 
  count(category, name = "n_trait_pairs")

write_csv(trait_pair_classification, "results/tab/trait_pair_classification_leaf_roots.csv")


# 10. Plotting function ----
plot_trait_relationship <- function(data, trait_pair, sig_general, sig_species, 
                                    species_colors, excluded_species = NULL) {
  
  # Filter data for this trait pair
  plot_data <- data |> filter(traits == trait_pair)
  
  # Check if general relationship is significant
  is_sig_general <- trait_pair %in% sig_general$traits
  
  # Get significant species for this trait pair
  sig_spp <- sig_species |> 
    filter(traits == trait_pair) |> 
    pull(species)
  
  # Extract trait names for axis labels
  trait_names <- str_split(trait_pair, "~")[[1]]
  y_var <- trait_names[1]
  x_var <- trait_names[2]
  
  # Create nice labels (customize as needed)
  label_lookup <- c(
    "srl" = "SRL (m g^-1)",
    "rd" = "RD (mm)",
    "rtd" = "RTD (g cm^-3)",
    "rdmc" = "RDMC (mg g^-1)",
    "bi" = "BI",
    "root_depth" = "Root depth (cm)",
    "sla" = "SLA (cm² g^-1)",
    "ldmc" = "LDMC (mg g^-1)",
    "leaf_thickness" = "Leaf thickness (mm)",
    "bgb_agb" = "BGB:AGB",
    "veg_height" = "Vegetation height (cm)"
  )
  
  x_label <- label_lookup[x_var]
  y_label <- label_lookup[y_var]
  
  # Determine plot subtitle based on significance
  subtitle_text <- case_when(
    is_sig_general & length(sig_spp) > 0 ~ 
      paste0("Significant: Overall + ", length(sig_spp), " species"),
    is_sig_general & length(sig_spp) == 0 ~ 
      "Significant: Overall only",
    !is_sig_general & length(sig_spp) > 0 ~ 
      paste0("Significant: ", length(sig_spp), " species only"),
    TRUE ~ "Not significant"
  )
  
  # Start building plot
  p <- ggplot(plot_data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_point(color = "grey", alpha = 0.6, size = 2)
  
  # Add general relationship line if significant
  if (is_sig_general) {
    p <- p + 
      geom_line(aes(y = predicted_all, color = "All plants"), 
                linewidth = 2, linetype = 1)
  }
  
  # Add species-specific lines if significant
  if (length(sig_spp) > 0) {
    species_data <- plot_data |> 
      filter(species %in% sig_spp)
    
    # Remove excluded species if specified
    if (!is.null(excluded_species)) {
      species_data <- species_data |> 
        filter(!species %in% excluded_species)
    }
    
    if (nrow(species_data) > 0) {
      p <- p +
        geom_line(data = species_data,
                  aes(y = predicted_spp, color = species, group = species),
                  linewidth = 1.5, linetype = 1)
    }
  }
  
  # Add color scale
  colors <- c("All plants" = "black", species_colors)
  p <- p + scale_color_manual(values = colors)
  
  # Add labels and theme
  p <- p +
    labs(x = x_label, y = y_label, 
         title = trait_pair,
         subtitle = subtitle_text) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 13, color = "gray30"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 12, color = "black"),
      axis.ticks = element_line(linewidth = 1),
      panel.border = element_blank(),
      axis.line = element_line(linewidth = 1, colour = "black"),
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      legend.position = "right"
    )
  
  return(p)
}


## Create all plots ----
all_plots <- map(
  all_sig_trait_pairs,
  ~ plot_trait_relationship(
    data = results_general,
    trait_pair = .x,
    sig_general = sig_best,
    sig_species = sig_best_spp,
    species_colors = species_colors,
    excluded_species = NULL
  )
)

# Name the list elements
names(all_plots) <- all_sig_trait_pairs

## Create separate plot lists by category ----
plots_both <- all_plots[trait_pair_classification |> 
                          filter(category == "Both") |> 
                          pull(traits)]

plots_general_only <- all_plots[trait_pair_classification |>
                                  filter(category == "General only") |>
                                  pull(traits)]

plots_species_only <- all_plots[trait_pair_classification |> 
                                  filter(category == "Species only") |> 
                                  pull(traits)]

## 11. Export plots -----

leaf_roots_pairwise_spp <- plot_grid(plots_species_only$`rd~leaf_thickness`, plots_species_only$`leaf_thickness~rd`,
                                     plots_species_only$`rtd~ldmc`, plots_species_only$`ldmc~rtd`,
                                     plots_species_only$`sla~srl`, plots_species_only$`srl~sla`,
                                     ncol = 2)

ggsave("results/img/pairwise/leaf_roots_pairwise_spp.png", leaf_roots_pairwise_spp, dpi = 300,
       width = 30, height = 30, units = "cm")
