################################################################################
# Trait-trait covariation supplementary data
################################################################################
#
# JDMW
# January, 2026
#
# Description: This code analyses trait-trait covariation among selected leaf,
# roots, and plant size traits
################################################################################

# SET UP #######################################################################

#  Load packages ----------------------------------------------------------------
pkgs <- c("lme4", "lmerTest", "broom.mixed", "modelsummary", "tidyverse", "glue", "fixest",
          "patchwork")
# lapply(pkgs, install.packages, character.only = TRUE)
lapply(pkgs, library, character.only = TRUE)
remove(pkgs)

# 1. Load data ----
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
traits <- c("belowground_biomass", "aboveground_biomass", "root_depth", "veg_height", "srl", "sla")

trait_data <- read_csv("data/processed/v_PFCT7_clean_functional_traits_2023.csv")
trait_data_wide <- read_csv("data/processed/v_PFCT7_clean_functional_traits_2023_wide.csv") %>%
  filter(aboveground_biomass < 40) %>% # filter out large AG measure
  mutate(across(all_of(traits), ~ as.numeric(range01(.x)))) # scale all vars
names(trait_data_wide)

# PROCESSING ###################################################################
# 3. Define traits ----
metadata <- c("id", "aspect", "site_id", "elevation_m_asl", "plant_id", "species", "family", "growth_form")

# 4. Generate ALL possible pairs----
pairs <- expand_grid(trait_x = traits, trait_y = traits) |> 
  filter(trait_x != trait_y) |>  # Remove self-pairs
  mutate(traits = paste(trait_y, trait_x, sep = "~")) |> 
  filter(traits %in% c(
    "aboveground_biomass~belowground_biomass", "veg_height~root_depth", 'root_depth~belowground_biomass','srl~belowground_biomass', 'sla~belowground_biomass', 'veg_height~aboveground_biomass','srl~aboveground_biomass', 'sla~aboveground_biomass'))

# Define grouping factor
group <- "species"  

# Express formulas
formulas <- tribble(
  ~type, ~formula_template,
  "linear", "{trait_y} ~ {trait_x} + (1 + {trait_x} | {group})",
  "poly", "{trait_y} ~ {trait_x} + I({trait_x}^2) + (1 + {trait_x} | {group})",
  "exp", "log({trait_y} + 1e-10) ~ {trait_x} + (1 + {trait_x} | {group})"
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
             REML = FALSE)),
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
                         predict(.x)
                       } else if (.y == "exp") {
                         exp(predict(.x))
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
    filter(traits == trait_pair)
  
  global_data <- global_plot_data |> 
    filter(traits == trait_pair)
  
  plot <- ggplot(spp_data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    ## Raw data
    geom_point(aes(color = species),alpha = 0.6,size = 2) +
    ## Species-level predictions
    geom_line(aes(y = predicted,color = species,linetype = significant),linewidth = 1) +
    ## Global prediction
    geom_line(data = global_data,aes(x = .data[[x_var]],y = predicted,linetype = significant), linewidth = 1.2,color = "black") +
    scale_color_manual(values = c("#41049EFF","#8104A7FF","#B6308BFF","#F79143FF","#FCCE15FF"),
                       name = "Species") +
    scale_linetype_manual(values = c("TRUE" = 1, "FALSE" = 2), guide = 'none') +
    labs(x = x_label, y = y_label) +
    theme_classic()
  
  result_plots[[i]] <- plot
}

# make a combined plot of all plots
PS_legend <- cowplot::get_legend(result_plots[[1]])

ps_combo_plot <- result_plots[[1]] +
  geom_abline(slope = 1, linewidth = 1.1, col = 'red') + result_plots[[8]] + cowplot::ggdraw(PS_legend) + result_plots[[2]] + result_plots[[3]] + result_plots[[4]] + 
  result_plots[[5]] + result_plots[[6]] + result_plots[[7]] +
  # plot_layout(guides = 'collect') &
  plot_annotation(tag_levels = list(c('A.','B.','', 'C.', 'D.', 'E.', 'F.', 'G.','H.'))) &
  theme(legend.position = 'none') 
ps_combo_plot

ggsave(filename = 'mixed_models_trait_trait_PS.png',
       path = 'results/img/pairwise/',
       width = 9.59, height = 8.00, dpi = 310)

#### Tidy up model summaries ----
#### a) Global model summaries ----
model_outputs <- modelsummary(results_best$model, 
                              estimate = "{estimate} [{conf.low}, {conf.high}] {stars}",
                              output = 'data.frame')
names(model_outputs)[4:11] <- results_best$model_name

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
names(model_outputs_clean) <- str_replace(names(model_outputs_clean), '_', ' ')
names(model_outputs_clean) <- str_replace(names(model_outputs_clean), 'poly', '(polynomial)')
names(model_outputs_clean) <- str_replace(names(model_outputs_clean), 'linear', '(linear)')
names(model_outputs_clean) <- str_replace(names(model_outputs_clean), 'exp', '(exponential)')

# change column orders, drop part
model_outputs_clean <- model_outputs_clean %>%
  select(term, 7, 11, 8, 10, 9, 6, 5, 4)

# Add letters on that match to the figure
names(model_outputs_clean)[2:9] <- paste0(LETTERS[1:8], '. ', names(model_outputs_clean)[2:9])

# Transpose format
model_outputs_wide <- model_outputs_clean %>%
  pivot_longer(cols = 2:9, names_to = 'model') %>%
  pivot_wider(names_from = 'term', values_from = 'value')

writexl::write_xlsx(model_outputs_wide, 'results/tab/trait_trait_PS_mixed_model_outputs_raw_GLOBAL.xlsx')

#### b) Species model summaries ----
spp_model_outputs <- modelsummary(results_best_spp$model, 
                                  estimate = "{estimate} [{conf.low}, {conf.high}] {stars}",
                                  output = 'data.frame')
names(spp_model_outputs)[4:43] <- paste(results_best_spp$model_name, results_best_spp$species)

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
names(spp_model_outputs_clean) <- str_replace(names(spp_model_outputs_clean), '_', ' ')
names(spp_model_outputs_clean) <- str_replace(names(spp_model_outputs_clean), 'poly', '(polynomial)')
names(spp_model_outputs_clean) <- str_replace(names(spp_model_outputs_clean), 'linear', '(linear)')
names(spp_model_outputs_clean) <- str_replace(names(spp_model_outputs_clean), 'exp', '(exponential)')

# change column orders, drop part
names(spp_model_outputs_clean)

spp_model_outputs_clean <- spp_model_outputs_clean %>%
  select(term, 
         3+4, 3+8+4, 3+8+8+4, 3+8+8+8+4, 3+8+8+8+8+4, # agb~bgb
         3+8, 3+8+8, 3+8+8+8, 3+8+8+8+8, 3+8+8+8+8+8, # vheight~rdepth
         3+5, 3+8+5, 3+8+8+5, 3+8+8+8+5, 3+8+8+8+8+5, # rdepth~bgb
         3+7, 3+8+7, 3+8+8+7, 3+8+8+8+7, 3+8+8+8+8+7, # srl~bgb
         3+6, 3+8+6, 3+8+8+6, 3+8+8+8+6, 3+8+8+8+8+6, # sla~bgb
         3+3, 3+8+3, 3+8+8+3, 3+8+8+8+3, 3+8+8+8+8+3, # vheight~agb
         3+2, 3+8+2, 3+8+8+2, 3+8+8+8+2, 3+8+8+8+8+2, # srl~agb
         3+1, 3+8+1, 3+8+8+1, 3+8+8+8+1, 3+8+8+8+8+1 # sla~agb
  )

# Add letters on that match to the figure
names(spp_model_outputs_clean)[2:41] <- paste0(rep(LETTERS[1:8], each = 5), '. ', names(spp_model_outputs_clean)[2:41])

# Transpose format
spp_model_outputs_wide <- spp_model_outputs_clean %>%
  pivot_longer(cols = 2:41, names_to = 'model') %>%
  pivot_wider(names_from = 'term', values_from = 'value') %>%
  mutate(species = str_extract(model, "\\S+\\s+\\S+$"),
         model = str_remove(model, "\\s*\\S+\\s+\\S+$")) %>%
  dplyr::select(1, 11, 1:10)

writexl::write_xlsx(spp_model_outputs_wide, 'results/tab/trait_trait_PS_mixed_model_outputs_raw_SPECIES.xlsx')
