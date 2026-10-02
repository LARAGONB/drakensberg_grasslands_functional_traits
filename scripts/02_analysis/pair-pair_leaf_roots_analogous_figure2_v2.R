################################################################################
# Trait-trait covariation
################################################################################
#
# JDMW
# January, 2026
#
# Description: This code analyses trait-trait covariation among analogous leaf
# and roots traits using linear, polynomial and exponential models
################################################################################

# SET UP #######################################################################

#  1. Load packages ----------------------------------------------------------------
pkgs <- c("lme4", "lmerTest", "broom.mixed", "modelsummary", "tidyverse", "glue", "fixest",
          "patchwork")
# lapply(pkgs, install.packages, character.only = TRUE)
lapply(pkgs, library, character.only = TRUE)
remove(pkgs)

# 2. Load data ----
range01 <- function(x) {(x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))}
traits <- c("leaf_thickness", "sla", "ldmc", "rd", "srl", "rtd", "belowground_biomass", "aboveground_biomass", "root_depth", "veg_height", "reproductive_height")

trait_data <- read_csv("data/processed/v_PFCT7_clean_functional_traits_2023.csv")


# PROCESSING ###################################################################
trait_data_wide <- read_csv("data/processed/v_PFCT7_clean_functional_traits_2023_wide.csv") %>%
  mutate(aboveground_biomass = case_when(
    aboveground_biomass > 40 ~ NA,
    TRUE ~ aboveground_biomass)) %>% # filter out large AG measure
  mutate(across(all_of(traits), ~ as.numeric(range01(.x)))) # scale all vars
  
# 3. Define traits ----
metadata <- c("id", "aspect", "site_id", "elevation_m_asl", "plant_id", "species", "family", "growth_form")

# 4. Generate ALL possible pairs----
pairs <- expand_grid(trait_x = traits, trait_y = traits) |> 
  filter(trait_x != trait_y) |>  # Remove self-pairs
  mutate(traits = paste(trait_y, trait_x, sep = "~")) |> 
  filter(traits %in% c("leaf_thickness~sla", "ldmc~sla", 'ldmc~leaf_thickness',
                       'rd~srl', 'rtd~srl', 'rtd~rd',
                       'sla~srl', 'leaf_thickness~rd', 'ldmc~rtd',
                       "aboveground_biomass~belowground_biomass", "veg_height~root_depth", "reproductive_height~veg_height"))

# Define grouping factor
group <- "species"  

# Express formulas
formulas <- tribble(
  ~type, ~formula_template,
  "linear", "{trait_y} ~ {trait_x} + (1 + {trait_x} | {group})",
  "poly", "{trait_y} ~ {trait_x} + I({trait_x}^2) + (1 + {trait_x} | {group})",
  "exp", "sign({trait_y}) * log(abs({trait_y}) + 1e-10) ~ {trait_x} + (1 + {trait_x} | {group})"
  # "exp", "log({trait_y} + 1e-10) ~ {trait_x} + (1 + {trait_x} | {group})"
)

# POPULATION/GLOBAL LEVEL ----
# Create the model grid combinations
model_grid <- crossing(pairs, formulas) |> 
  mutate(
    base_formula = pmap_chr(
      list(formula_template, trait_x, trait_y, group),
      ~ glue(
        ..1,
        trait_x = ..2,
        trait_y = ..3,
        group   = ..4
      )
    ),
    model_name = paste(traits, type, sep = "_")
  ) |> 
  select(!formula_template)

# Run GLOBAL models
results <- model_grid |> 
  mutate(
    model = map(
      base_formula,
      ~ lmer(as.formula(.x),
             data = trait_data_wide,
             REML = FALSE,  
             na.action = na.exclude,
             control = lmerControl(
               optimizer = "bobyqa",
               optCtrl = list(maxfun = 2e5),
               check.conv.singular = "ignore"
             )
             )),
    p_value = map2_dbl(model, trait_x, ~ summary(.x)$coefficients[.y, "Pr(>|t|)"]),
    significant = p_value < 0.05,
    aic  = map_dbl(model, AIC),
    predicted = map2(model, type, ~ if (.y %in% c("linear", "poly")) {predict(.x, re.form = NA)} else {
        exp(predict(.x, re.form = NA))}))

# Sort table based on aic value for each trait pair ----
results_sort <- results |> 
  group_by(traits) |> 
  arrange(traits,aic)

# Extract best model (lowest aic) for each trait pair ----
results_best <- results %>%
  group_by(traits) %>%
  filter(aic == min(aic)) %>%
  ungroup()

# SPECIES LEVEL ----
## Nest data by species
nested_spp <- trait_data_wide |> 
  group_by(species) |> 
  nest()

## Expand grid to species × model grid
species_grid <- crossing(
  species = unique(trait_data_wide$species),
  model_grid
) |> 
  left_join(nested_spp, by = "species")

### run SPECIES models
results_spp <- species_grid |> 
  mutate(
    ## Strip random effects from formula
    fixed_formula = str_remove(base_formula, "\\+ \\(.*\\)"),
    model = map2(fixed_formula, data, ~ feols(as.formula(.x), data = .y)),
    p_value = map_dbl(model, ~ fitstat(.x, "f")$f$p, se = "hc3"),
    significant = p_value < 0.05,
    aic = map_dbl(model, ~ fitstat(.x, "aic")$aic, se = "hc3"),
    predicted = map2(model, type, ~ 
                       if (.y %in% c("linear", "poly")) {
                         predict(.x, sample = 'original')
                       } else if (.y == "exp") {
                         exp(predict(.x, sample = 'original'))
                       } else {
                         NA  # or handle other model types
                       }))

results_best_spp <- results_spp |> 
  group_by(species, traits) |> 
  filter(aic == min(aic)) |> 
  ungroup()

spp_plot_data <- results_best_spp |> 
  mutate(
    data = map2(data, predicted, ~ mutate(.x, predicted = .y))
  ) |> 
  select(species, traits, significant, data) |> 
  unnest(data)

global_plot_data <- results_best |> 
  mutate(
    data = map(predicted, ~ mutate(trait_data_wide, predicted = .x))
  ) |> 
  select(traits, significant, data) |> 
  unnest(data)

#Species order
spp_plot_data$species <- factor(spp_plot_data$species, levels = c("Eragrostis capensis", "Harpochloa falx", "Themeda triandra", "Helichrysum pilosellum",  "Senecio glaberrimus")) 
spp_plot_data <- spp_plot_data %>%
  mutate(species = case_when(
    species == "Eragrostis capensis" ~ "ERCA",
    species == "Harpochloa falx" ~ "HAFA",
    species == "Themeda triandra" ~ "THTR",
    species == "Helichrysum pilosellum" ~ "HEPI",
    species == "Senecio glaberrimus" ~ "SEGL"
  ))

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
  "veg_height" = "Vegetation height (cm)",
  "reproductive_height" = "Reproductive height (cm)",
  "belowground_biomass" = 'Belowground biomass (g)',
  "aboveground_biomass" = 'Aboveground biomass (g)'
)

result_plots <- list()

for(i in 1:length(pairs$traits)){
  trait_pair <- pairs$traits[[i]]
  x_var <- str_split(trait_pair, '~')[[1]][2]
  y_var <- str_split(trait_pair, '~')[[1]][1]

  x_label <- label_lookup[x_var]
  y_label <- label_lookup[y_var]
  
  spp_data <- spp_plot_data |> 
    filter(traits == trait_pair) |> 
    filter(trait_pair != "aboveground_biomass~belowground_biomass" |
             (!is.na(aboveground_biomass) & !is.na(belowground_biomass))) |> 
    dplyr::mutate(
      species = factor(species, levels = c("ERCA", "HAFA", "THTR", "HEPI", "SEGL")))
  
  global_data <- global_plot_data |> 
    filter(traits == trait_pair) |> 
    filter(trait_pair != "aboveground_biomass~belowground_biomass" |
             (!is.na(aboveground_biomass) & !is.na(belowground_biomass)))
  
  plot <- ggplot(spp_data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_point(aes(color = species, shape = species), alpha = 0.6, size = 2) +
    geom_line(aes(y = predicted, color = species, linetype = significant), linewidth = 1) +
    geom_line(
      data = global_data,
      aes(x = .data[[x_var]], y = predicted, linetype = significant),
      linewidth = 1.2, color = "black"
    ) +
    scale_color_manual(
      values = c("#42049EFF","#8204A7FF","#B6308BFF","#F79143FF","#FCCE25FF"),
      name = "Species"
    ) +
    scale_shape_manual(
      values = c("ERCA" = 16, "HAFA" = 16, "THTR" = 16, "HEPI" = 15, "SEGL" = 15),
      name = "Species"
    ) +
    scale_linetype_manual(values = c("TRUE" = 1, "FALSE" = 2), guide = "none") +
    labs(x = x_label, y = y_label) +
    theme_classic()
  
  result_plots[[i]] <- plot
}

result_plots[[11]]

# make a combined plot of all plots
combo_plot <- result_plots[[2]] + result_plots[[3]] + result_plots[[1]] + 
  result_plots[[7]] + result_plots[[8]] + result_plots[[5]] +
  result_plots[[6]] + result_plots[[4]] + result_plots[[9]] +
  plot_layout(guides = 'collect') &
  plot_annotation(tag_levels = 'A', tag_suffix = '.')
combo_plot

ggsave(filename = 'mixed_models_trait_trait.png',
       path = 'results/img/pairwise/',
       width = 9.59, height = 8.00, dpi = 320)

# make a combined plot of all plots
analag_plot <- 
  result_plots[[6]] + result_plots[[4]] + result_plots[[9]] +
  result_plots[[10]]+
  geom_abline(slope = 1, linewidth = 1.1, col = 'red') +
  
  plot_layout(guides = 'collect') &
  plot_annotation(tag_levels = 'A', tag_suffix = '.')
analag_plot

ggsave(filename = 'mixed_models_trait_trait_analog.png',
       path = 'results/img/pairwise/',
       width = 9.59, height = 2.8, dpi = 320)


#### updated plot with AGB-BGB and VH and RD
PS_legend <- cowplot::get_legend(result_plots[[1]])

# make a combined plot of all plots
analag_plot <- 
  result_plots[[6]] + result_plots[[4]] + result_plots[[9]] +
  result_plots[[10]]+
  geom_abline(slope = 1, linewidth = 1.1, col = 'red') +
  result_plots[[11]] + result_plots[[12]] +
  # cowplot::ggdraw(PS_legend) + 
  plot_annotation(tag_levels = 'A', tag_suffix = '.') +
  plot_layout(guides = 'collect') 
analag_plot

ggsave(filename = 'mixed_models_trait_trait_analog_v2.png',
       path = 'results/img/pairwise/',
       width = 9.59, height = 5.6, dpi = 320)


# make a combined plot of all plots
leaf_root_plot <- 
  result_plots[[2]] + result_plots[[3]] + result_plots[[1]] + 
  result_plots[[7]] + result_plots[[8]] + result_plots[[5]] +
  plot_layout(guides = 'collect') &
  plot_annotation(tag_levels = 'A', tag_suffix = '.')
leaf_root_plot

ggsave(filename = 'mixed_models_trait_trait_leaf_root_v2.png',
       path = 'results/img/pairwise/',
       width = 9.59, height = 2.8*2, dpi = 320)

#### Tidy up model summaries ----
#### a) Global model summaries ----
model_outputs <- modelsummary(results_best$model, 
             estimate = "{estimate} [{conf.low}, {conf.high}] {stars}",
             output = 'data.frame')
names(model_outputs)[4:15] <- results_best$model_name

model_outputs_clean <- model_outputs %>%
  mutate(term = str_replace_all(term, paste(traits, collapse = '|'), 'trait_x')) %>% # replace trait names with trait_x
  pivot_longer(cols = !c('part', 'term', 'statistic')) %>%
  mutate(statistic = case_when(
    part == 'gof' ~ 'statistic',
    statistic == 'std.error' ~ 'SD',
    TRUE ~ statistic
  ),
  value = case_when(
    value == '' ~ NA_character_,
    TRUE ~ value
  )) %>%
  filter(!is.na(value)) %>% # remove empty cells
  pivot_wider(names_from = 'name', values_from = 'value') %>%
  mutate(part = case_when(
    part == 'gof' ~ 'goodness-of-fit',
    TRUE ~ part
  ),
  term = case_when(
    term == '(Intercept)' ~ 'Intercept',
    term == 'I(trait_x^2)' ~ 'trait_x^2',
    TRUE ~ term
  )) %>%
  filter(statistic != 'SD') # filter out SD as we have CI

# remove underscores
names(model_outputs_clean) <- str_replace(names(model_outputs_clean), '_', ' ')
names(model_outputs_clean) <- str_replace(names(model_outputs_clean), '_', ' ')
names(model_outputs_clean) <- str_replace(names(model_outputs_clean), 'poly', '(polynomial)')
names(model_outputs_clean) <- str_replace(names(model_outputs_clean), 'linear', '(linear)')
names(model_outputs_clean) <- str_replace(names(model_outputs_clean), 'exp', '(exponential)')

# change column orders, drop part
model_outputs_clean <- model_outputs_clean %>%
  select(term, 9, 8, 4, 10, 11, 6, 12, 5, 7)

# Add letters on that match to the figure
names(model_outputs_clean)[2:10] <- paste0(LETTERS[1:9], '. ', names(model_outputs_clean)[2:10])

# Transpose format
model_outputs_wide <- model_outputs_clean %>%
  pivot_longer(cols = 2:10, names_to = 'model') %>%
  pivot_wider(names_from = 'term', values_from = 'value')

writexl::write_xlsx(model_outputs_wide, 'results/tab/trait_trait_mixed_model_outputs_raw_GLOBAL.xlsx')

#### b) Species model summaries ----
spp_model_outputs <- modelsummary(results_best_spp$model, 
                              estimate = "{estimate} [{conf.low}, {conf.high}] {stars}",
                              output = 'data.frame')
names(spp_model_outputs)[4:48] <- paste(results_best_spp$model_name, results_best_spp$species)

spp_model_outputs_clean <- spp_model_outputs %>%
  mutate(term = str_replace_all(term, paste(traits, collapse = '|'), 'trait_x')) %>% # replace trait names with trait_x
  pivot_longer(cols = !c('part', 'term', 'statistic')) %>%
  mutate(statistic = case_when(
    part == 'gof' ~ 'statistic',
    statistic == 'std.error' ~ 'SD',
    TRUE ~ statistic
  ),
  value = case_when(
    value == '' ~ NA_character_,
    TRUE ~ value
  )) %>%
  filter(!is.na(value)) %>% # remove empty cells
  pivot_wider(names_from = 'name', values_from = 'value') %>%
  mutate(part = case_when(
    part == 'gof' ~ 'goodness-of-fit',
    TRUE ~ part
  ),
  term = case_when(
    term == '(Intercept)' ~ 'Intercept',
    term == 'I(I(trait_x^2))' ~ 'trait_x^2',
    TRUE ~ term
  )) %>%
  filter(statistic != 'SD') %>%
  filter(term != 'Std.Errors')

# remove underscores
names(spp_model_outputs_clean) <- str_replace(names(spp_model_outputs_clean), '_', ' ')
names(spp_model_outputs_clean) <- str_replace(names(spp_model_outputs_clean), '_', ' ')
names(spp_model_outputs_clean) <- str_replace(names(spp_model_outputs_clean), 'poly', '(polynomial)')
names(spp_model_outputs_clean) <- str_replace(names(spp_model_outputs_clean), 'linear', '(linear)')
names(spp_model_outputs_clean) <- str_replace(names(spp_model_outputs_clean), 'exp', '(exponential)')

# change column orders, drop part
names(spp_model_outputs_clean)

spp_model_outputs_clean <- spp_model_outputs_clean %>%
  select(term, 
         3+6, 3+9+6, 3+9+9+6, 3+9+9+9+6, 3+9+9+9+9+6, # lt~sla
         3+5, 3+9+5, 3+9+9+5, 3+9+9+9+5, 3+9+9+9+9+5, # ldmc~sla
         3+1, 3+9+1, 3+9+9+1, 3+9+9+9+1, 3+9+9+9+9+1, # ldmc~lt
         3+7, 3+9+7, 3+9+9+7, 3+9+9+9+7, 3+9+9+9+9+7, # rd~srl
         3+8, 3+9+8, 3+9+9+8, 3+9+9+9+8, 3+9+9+9+9+8, # rtd~srl
         3+3, 3+9+3, 3+9+9+3, 3+9+9+9+3, 3+9+9+9+9+3, # rtd~rd
         3+9, 3+9+9, 3+9+9+9, 3+9+9+9+9, 3+9+9+9+9+9, # sla~srl
         3+2, 3+9+2, 3+9+9+2, 3+9+9+9+2, 3+9+9+9+9+2, # lt~rd
         3+4, 3+9+4, 3+9+9+4, 3+9+9+9+4, 3+9+9+9+9+4, # ldmc~rtd
  )

# Add letters on that match to the figure
names(spp_model_outputs_clean)[2:46] <- paste0(rep(LETTERS[1:9], each = 5), '. ', names(spp_model_outputs_clean)[2:46])

# Transpose format
spp_model_outputs_wide <- spp_model_outputs_clean %>%
  pivot_longer(cols = 2:46, names_to = 'model') %>%
  pivot_wider(names_from = 'term', values_from = 'value') %>%
  mutate(species = str_extract(model, "\\S+\\s+\\S+$"),
    model = str_remove(model, "\\s*\\S+\\s+\\S+$")) %>%
  dplyr::select(1, 11, 2:10)

writexl::write_xlsx(spp_model_outputs_wide, 'results/tab/trait_trait_mixed_model_outputs_raw_SPECIES.xlsx')