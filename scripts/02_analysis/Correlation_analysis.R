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

#install.packages("devtools", "tidyverse", "ggplot2", "plotly", "vegan", "ggvegan", "ggrepel", "glue", "viridis", "fixest", "lmtest", "corrplot", "skedastic") #install if needed
# devtools::install_github("gavinsimpson/ggvegan")
pkgs <- c("devtools", "tidyverse", "ggplot2", "plotly", "vegan", "ggvegan", "ggrepel", "glue", "viridis", "fixest",
          "lmtest", "corrplot", "skedastic", "patchwork")
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
pval_matrix_traits <- cor.mtest(matrix_traits, use = "pairwise.complete.obs", method = "pearson")

## Visual representation ----
png("results/img/correlation_plot.png", width = 15, heigh = 15, unit = "cm", res = 300)
corrplot(
  cor_matrix_traits,
  method = "color",
  type = "lower",
  tl.srt = 45,
  addCoef.col = "black",      # Show coefficients
  tl.col = "black",           # Axis text color
  p.mat = pval_matrix_traits$p, # P-value matrix
  sig.level = 0.05,           # Significance threshold
  insig = "blank",
  pch.col = "white",
  pch.cex = 4,
  number.digits = 2           # Number of decimals for coefficients
)
dev.off()
# corrplot(
#   cor_matrix_traits,
#   method = "color",
#   type = "lower",
#   tl.srt = 45,
#   addCoef.col = "black",      # Show coefficients
#   tl.col = "black",           # Axis text color
#   p.mat = pval_matrix_traits$p, # P-value matrix
#   sig.level = 0.05,           # Significance threshold
#   insig = "pch",
#   pch = 4,
#   pch.col = "white",
#   pch.cex = 4,
#   number.digits = 2           # Number of decimals for coefficients
# )

# 4. Correlations across species using linear, polynomial and exponential model ----

## Trait pairs ----
pairs <- tribble(
  ~trait_x, ~trait_y,
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
  "bgb_agb", "ldmc") |> 
  mutate(traits = paste(trait_y,trait_x,sep = "~"))


## Model types ----
formulas <- tribble(
  ~type, ~formula_template,
  "linear", "{trait_y} ~ {trait_x}",
  "poly", "{trait_y} ~ {trait_x} + I({trait_x}^2)",
  "exp", "log({trait_y}) ~ {trait_x}"
)

## Expand for all pairs and formulas ----
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

## Fit all models and extract valuable information ----
results <- model_grid |> 
  mutate(
    model = map(base_formula, ~ feols(as.formula(.x), data = trait_data_wide)),
    estimate = map_dbl(model, ~ coef(.x)[[2]]),
    mean_dv = map_dbl(model, ~ fitstat(.x, "my")$my),
    # se = map_dbl(model, ~ fixest::se(.x)[[2]]),
    rmse = map_dbl(model, ~ fitstat(.x, "rmse")$rmse),
    f_value = map_dbl(model, ~ fitstat(.x, "f")$f$stat),
    p_value = map_dbl(model, ~ fitstat(.x, "f")$f$p),
    ar2 = map_dbl(model, ~ fitstat(.x, "ar2")$ar2),
    significant = p_value < 0.05,
    aic = map_dbl(model, ~ fitstat(.x, "aic")$aic),
    w_shapiro = map_dbl(model, ~ shapiro.test(resid(.x))$statistic[[1]]),
    p_shapiro = map_dbl(model, ~ shapiro.test(resid(.x))$p.value),
    w_shap_normal = w_shapiro > 0.9,
    p_shap_normal = p_shapiro >= 0.05,
    s_white = map_dbl(base_formula, ~ white(lm(as.formula(.x), data = trait_data_wide))$statistic),
    p_white = map_dbl(base_formula, ~ white(lm(as.formula(.x), data = trait_data_wide))$p.value),
    p_whit_normal = p_white >= 0.05,
    p_bp = map_dbl(model, ~bptest(resid(.x) ~ fitted(.x))$p.value),
    p_bp_normal = p_bp >= 0.05,
    etable = map(model, ~ etable(.x, fitstat = ~ . + f + my + rmse)),
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

sig_best

## Extract the fitted values per trait pair ----
selected_cols <- c(
  "id", "aspect", "site_id", "elevation_m_asl", "plant_id", "species", "family", "growth_form",
  "rd", "root_depth", "bi", "rtd", "rdmc", "ldmc", "bgb_agb", "srl", "sla", "predicted"
)

results_blong <- results_best |> 
  mutate(
    data = map(predicted, ~ mutate(trait_data_wide, predicted = .x))
  )  |> 
  select(1:6, data) |> 
  unnest(data) |> 
  select(1:6, all_of(selected_cols)) |> 
  rename(predicted_all = predicted)

# 5. Correlations within species ----

# Nest data by species
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
    estimate = map_dbl(model, ~ coef(.x)[[2]]),
    mean_dv = map_dbl(model, ~ fitstat(.x, "my")$my),
    se = map_dbl(model, ~ fixest::se(.x)[[2]]),
    rmse = map_dbl(model, ~ fitstat(.x, "rmse")$rmse),
    f_value = map_dbl(model, ~ fitstat(.x, "f")$f$stat),
    p_value = map_dbl(model, ~ fitstat(.x, "f")$f$p),
    ar2 = map_dbl(model, ~ fitstat(.x, "ar2")$ar2),
    significant = p_value < 0.05,
    aic = map_dbl(model, ~ fitstat(.x, "aic")$aic),
    w_shapiro = map_dbl(model, ~ shapiro.test(resid(.x))$statistic[[1]]),
    p_shapiro = map_dbl(model, ~ shapiro.test(resid(.x))$p.value),
    w_shap_normal = w_shapiro > 0.9,
    p_shap_normal = p_shapiro >= 0.05,
    s_white = map_dbl(base_formula, ~ white(lm(as.formula(.x), data = trait_data_wide))$statistic),
    p_white = map_dbl(base_formula, ~ white(lm(as.formula(.x), data = trait_data_wide))$p.value),
    p_whit_normal = p_white >= 0.05,
    p_bp = map_dbl(model, ~bptest(resid(.x) ~ fitted(.x))$p.value),
    p_bp_normal = p_bp >= 0.05,
    etable = map(model, ~ etable(.x, fitstat = ~ . + f + my + rmse)),
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

sig_best_spp

## Extracte the fitted values per species and trait pair ----
selected_cols_spp <- c(
  "id", "aspect", "site_id", "elevation_m_asl", "plant_id", "family", "growth_form",
  "rd", "root_depth", "bi", "rtd", "rdmc", "ldmc", "bgb_agb", "srl", "sla", "predicted"
)

results_blong_spp <- results_best_spp %>%
  mutate(
    data = map2(data, predicted, ~ mutate(.x, predicted = .y))
  ) %>%
  select(1:7, data) %>%   # First 7 cols + data
  unnest(data) |> 
  select(1:7, all_of(selected_cols_spp)) |> 
  rename(predicted_spp = predicted)

# 7. Create full results table ----

join_cols <- setdiff(intersect(names(results_blong), names(results_blong_spp)), c("type", "base_formula", "model_name"))

results_general <- left_join(results_blong,results_blong_spp, by = join_cols) |> 
  mutate(species = fct_relevel(species, "Themeda triandra", after = 2))

# 6. Visualize ----
##Unique species colors

species_colors <- c(
  "Eragrostis capensis" = "#42049EFF",
  "Harpochloa falx" = "#8204A7FF",
  "Themeda triandra" = "#B6308BFF",
  "Helichrysum pilosellum" = "#F79143FF",
  "Senecio glaberrimus" = "#FCCE25FF")

## Using facet_wrap ----
plot_df <- map2_dfr(pairs$trait_x, pairs$trait_y, ~ {
  trait_data_wide |> 
    select(x = all_of(.x), y = all_of(.y)) |> 
    mutate(
      traits = paste(.y, .x, sep = "~"),
      trait_x = paste(.x),
      trait_y = paste(.y)) |> 
    relocate(x, .after = last_col()) |> 
    relocate(y, .after = last_col())
  }
  )

ggplot(plot_df, aes(x = x, y = y)) +
  geom_point() +
  facet_wrap(~ traits, scales = "free") +
  theme_bw()

## Unique plots and wrap (prefered) ----
#### srl ~ rtd ----

srl_rtd <- results_general |> 
  filter(traits == "srl~rtd") |> 
  ggplot(aes(x = rtd, y = srl)) +
  geom_point(color = "grey") +
  geom_line(aes(y = predicted_all, color = "All plants"), linewidth = 2, linetype = 1) +
  geom_line(aes(y = predicted_spp, color = species, group = species), linewidth = 1.5,
            linetype = 1) +
  scale_color_manual(values = c("All plants" = "black", species_colors)) +
  labs(
    x = expression("RTD (g cm"^-3*")"),
    y = expression("SRL (m g"^-1*")")) + 
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.ticks = element_line(linewidth = 1),
    panel.border = element_blank(), 
    axis.line = element_line(linewidth = 1, colour = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 16))

srl_rtd

#### srl ~ rd ----

srl_rd <- results_general |> 
  filter(traits == "srl~rd") |> 
  ggplot(aes(x = rd, y = srl)) +
  geom_point(color = "grey") +
  geom_line(aes(y = predicted_all, color = "All plants"), linewidth = 2, linetype = 1) +
  geom_line(data = results_general |> 
              filter(traits == "srl~rd",
                     !species %in% "Senecio glaberrimus"),
            aes(y = predicted_spp, color = species, group = species), linewidth = 1.5,
            linetype = 1) +
  scale_color_manual(values = c("All plants" = "black", species_colors)) +
  labs(
    x = expression("RD (mm)"),
    y = expression("SRL (m g"^-1*")")) + 
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.ticks = element_line(linewidth = 1),
    panel.border = element_blank(), 
    axis.line = element_line(linewidth = 1, colour = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 16))

srl_rd

#### sla ~ root_depth ----

sla_root_depth <- results_general |> 
  filter(traits == "sla~root_depth") |> 
  ggplot(aes(x = root_depth, y = sla)) +
  geom_point(color = "grey") +
  geom_line(aes(y = predicted_all, color = "All plants"), linewidth = 2, linetype = 1) +
  # geom_line(aes(y = predicted_spp, color = species, group = species), linewidth = 1.5,
  #           linetype = 1) +
  scale_color_manual(values = c("All plants" = "black")) +
  labs(
    x = expression("RDepth (cm)"),
    y = expression("SLA (cm"^2*" g"^-1*")")) + 
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.ticks = element_line(linewidth = 1),
    panel.border = element_blank(), 
    axis.line = element_line(linewidth = 1, colour = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 16))


sla_root_depth

#### sla ~ ldmc ----

sla_ldmc <- results_general |> 
  filter(traits == "sla~ldmc") |> 
  ggplot(aes(x = ldmc, y = sla)) +
  geom_point(color = "grey") +
  geom_line(aes(y = predicted_all, color = "All plants"), linewidth = 2, linetype = 1) +
  geom_line(data = results_general |> 
              filter(traits == "sla~ldmc",
                     !species %in% c("Themeda triandra","Helichrysum pilosellum")),
            aes(y = predicted_spp, color = species, group = species), linewidth = 1.5,
            linetype = 1) +
  scale_color_manual(values = c("All plants" = "black", species_colors)) +
  labs(
    x = expression("LDMC (g g"^-1*")"),
    y = expression("SLA (cm"^2*" g"^-1*")")) + 
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.ticks = element_line(linewidth = 1),
    panel.border = element_blank(), 
    axis.line = element_line(linewidth = 1, colour = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 16))

sla_ldmc

#### sla ~ bi ----

sla_bi <- results_general |> 
  filter(traits == "sla~bi") |> 
  ggplot(aes(x = bi, y = sla)) +
  geom_point(color = "grey") +
  # geom_line(aes(y = predicted_all, color = "All plants"), linewidth = 2, linetype = 1) +
  # geom_line(aes(y = predicted_spp, color = species, group = species), linewidth = 1.5,
  #           linetype = 1) +
  # scale_color_manual(values = c("All plants" = "black", species_colors)) +
  labs(
    x = expression("BI (count mm"^-1*")"),
    y = expression("SLA (cm"^2*" g"^-1*")")) + 
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.ticks = element_line(linewidth = 1),
    panel.border = element_blank(), 
    axis.line = element_line(linewidth = 1, colour = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 16))

sla_bi

#### rtd ~ root_depth ----

rtd_root_depth <- results_general |> 
  filter(traits == "rtd~root_depth") |> 
  ggplot(aes(x = root_depth, y = rtd)) +
  geom_point(color = "grey") +
  geom_line(aes(y = predicted_all, color = "All plants"), linewidth = 2, linetype = 1) +
  geom_line(data = results_general |> 
              filter(traits == "rtd~root_depth",
                     species %in% "Harpochloa falx"),
            aes(y = predicted_spp, color = species, group = species), linewidth = 1.5,
            linetype = 1) +
  scale_color_manual(values = c("All plants" = "black", species_colors)) +
  labs(
    x = expression("RDepth (cm)"),
    y = expression("RTD (g cm"^-3*")")) + 
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.ticks = element_line(linewidth = 1),
    panel.border = element_blank(), 
    axis.line = element_line(linewidth = 1, colour = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 16))

rtd_root_depth

#### rtd ~ rdmc ----

rtd_rdmc <- results_general |> 
  filter(traits == "rtd~rdmc") |> 
  ggplot(aes(x = rdmc, y = rtd)) +
  geom_point(color = "grey") +
  geom_line(aes(y = predicted_all, color = "All plants"), linewidth = 2, linetype = 1) +
  geom_line(data = results_general |> 
              filter(traits == "rtd~rdmc",
                     species %in% "Harpochloa falx"),
            aes(y = predicted_spp, color = species, group = species), linewidth = 1.5,
            linetype = 1) +
  scale_color_manual(values = c("All plants" = "black", species_colors)) +
  labs(
    x = expression("RDMC (mg g"^-1*")"),
    y = expression("RTD (g cm"^-3*")")) + 
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.ticks = element_line(linewidth = 1),
    panel.border = element_blank(), 
    axis.line = element_line(linewidth = 1, colour = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 16))

rtd_rdmc

#### rtd ~ bi ----

rtd_bi <- results_general |> 
  filter(traits == "rtd~bi") |> 
  ggplot(aes(x = bi, y = rtd)) +
  geom_point(color = "grey") +
  geom_line(aes(y = predicted_all, color = "All plants"), linewidth = 2, linetype = 1) +
  geom_line(data = results_general |> 
              filter(traits == "rtd~bi",
                     species %in% "Themeda triandra"),
            aes(y = predicted_spp, color = species, group = species), linewidth = 1.5,
            linetype = 1) +
  scale_color_manual(values = c("All plants" = "black", species_colors)) +
  labs(
    x = expression("BI (count mm"^-1*")"),
    y = expression("RTD (g cm"^-3*")")) + 
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.ticks = element_line(linewidth = 1),
    panel.border = element_blank(), 
    axis.line = element_line(linewidth = 1, colour = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 16))

rtd_bi

#### ldmc ~ root_depth ----

ldmc_root_depth <- results_general |> 
  filter(traits == "ldmc~root_depth") |> 
  ggplot(aes(x = root_depth, y = ldmc)) +
  geom_point(color = "grey") +
  geom_line(aes(y = predicted_all, color = "All plants"), linewidth = 2, linetype = 1) +
  # geom_line(aes(y = predicted_spp, color = species, group = species), linewidth = 1.5,
  #           linetype = 1) +
  scale_color_manual(values = c("All plants" = "black")) +
  labs(
    x = expression("RDepth (cm)"),
    y = expression("LDMC (g g"^-1*")")) + 
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.ticks = element_line(linewidth = 1),
    panel.border = element_blank(), 
    axis.line = element_line(linewidth = 1, colour = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 16))

ldmc_root_depth

#### ldmc ~ bi ----

ldmc_bi <- results_general |> 
  filter(traits == "ldmc~bi") |> 
  ggplot(aes(x = bi, y = ldmc)) +
  geom_point(color = "grey") +
  geom_line(aes(y = predicted_all, color = "All plants"), linewidth = 2, linetype = 1) +
  # geom_line(aes(y = predicted_spp, color = species, group = species), linewidth = 1.5,
  #           linetype = 1) +
  scale_color_manual(values = c("All plants" = "black")) +
  labs(
    x = expression("BI (count mm"^-1*")"),
    y = expression("LDMC (g g"^-1*")")) + 
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.ticks = element_line(linewidth = 1),
    panel.border = element_blank(), 
    axis.line = element_line(linewidth = 1, colour = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 16))


ldmc_bi

#### ldmc ~ bgb_agb ----

ldmc_bgb_agb <- results_general |> 
  filter(traits == "ldmc~bgb_agb") |> 
  ggplot(aes(x = bgb_agb, y = ldmc)) +
  geom_point(color = "grey") +
  geom_line(aes(y = predicted_all, color = "All plants"), linewidth = 2, linetype = 1) +
  geom_line(data = results_general |> 
              filter(traits == "ldmc~bgb_agb",
                     species %in% "Senecio glaberrimus"),
            aes(y = predicted_spp, color = species, group = species), linewidth = 1.5,
            linetype = 1) +
  scale_color_manual(values = c("All plants" = "black", species_colors)) +
  labs(
    x = expression("BG:AG (g g"^-1*")"),
    y = expression("LDMC (g g"^-1*")")) + 
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.ticks = element_line(linewidth = 1),
    panel.border = element_blank(), 
    axis.line = element_line(linewidth = 1, colour = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 16))

ldmc_bgb_agb

# 7. Combined plot ----

## Remove legends ----

ldmc_bgb_agb_n <- ldmc_bgb_agb + theme(legend.position = "none")
ldmc_bi_n <- ldmc_bi + theme(legend.position = "none")
ldmc_root_depth_n <- ldmc_root_depth + theme(legend.position = "none")
rtd_bi_n <- rtd_bi + theme(legend.position = "none")
rtd_rdmc_n <- rtd_rdmc + theme(legend.position = "none")
rtd_root_depth_n <- rtd_root_depth + theme(legend.position = "none")
sla_bi_n <- sla_bi + theme(legend.position = "none")
sla_ldmc_n <- sla_ldmc + theme(legend.position = "none")
sla_root_depth_n <- sla_root_depth + theme(legend.position = "none")
srl_rd_n <- srl_rd + theme(legend.position = "none")
srl_rtd_n <- srl_rtd + theme(legend.position = "none")
legend <- cowplot::get_legend(srl_rtd)

## Combine plot ----

combined_plot <- plot_grid(
  ldmc_bgb_agb_n, ldmc_bi_n, ldmc_root_depth_n,
  rtd_bi_n, rtd_rdmc_n, rtd_root_depth_n, sla_bi_n,
  sla_ldmc_n, sla_root_depth_n, srl_rd_n, srl_rtd_n, legend,
  ncol = 3)
combined_plot

## Save plots ----

ggsave("results/img/traits_regressions_plot_tiff.tiff", combined_plot,
       width = 30, height = 30, units = "cm", dpi = 300)
ggsave("results/img/traits_regressions_plot_png.png", combined_plot,
       width = 30, height = 30, units = "cm", dpi = 300)

