################################################################################
# Trait-trait covariation supplementary data
################################################################################
#
# JDMW
# January, 2026
#
# Description: This code analyses trait-trait covariation among leaf traits including
# data from the aboveground trait data set
################################################################################

# SET UP #######################################################################

# 1. Load packages ----------------------------------------------------------------
pkgs <- c("lme4", "lmerTest", "broom.mixed", "modelsummary", "tidyverse", "glue", "fixest",
          "patchwork")
# lapply(pkgs, install.packages, character.only = TRUE)
lapply(pkgs, library, character.only = TRUE)
remove(pkgs)

### 2. Load and format data ----
### retrieve raw data files from the OSF project page
# osf_retrieve_node('hk2cy') %>%
#   osf_ls_files(path = 'iv_aboveground_traits/') %>%
#   filter(!str_detect(name, 'experiment')) %>% # don't download the raw scan files
#   osf_download(path = 'data/raw/all_leaf_traits/', conflicts = 'overwrite')

### load in leaf trait field data
leaf <- read_csv('data/raw/all_leaf_traits/iv_PFTC7_clean_elevationgradient_traits_2023.csv')
names(leaf)
table(leaf$traits)
table(leaf$species)

### To remove outliers, find values that sit outside of the IQR:
# https://statsandr.com/blog/outliers-detection-in-r/
key_leaf_t <- leaf %>%
  filter(traits %in% c('ldmc', 'leaf_thickness', 'sla')) %>%
  dplyr::select(-unit) %>%
  pivot_wider(names_from = 'traits', values_from = 'value') %>%
  filter(is.na(problem_flag)) 

### filter to select traits 
clean_leaf <- leaf %>%
  filter(traits %in% c('ldmc', 'leaf_thickness', 'sla')) %>%
  dplyr::select(-unit) %>%
  pivot_wider(names_from = 'traits', values_from = 'value') %>%
  filter(is.na(problem_flag)) %>%
  drop_na(c(sla, leaf_thickness, ldmc)) # drop nas

### To remove outliers, find values that sit outside of the IQR:
# https://statsandr.com/blog/outliers-detection-in-r/
boxplot(clean_leaf$sla)
boxplot(clean_leaf$ldmc)
boxplot(clean_leaf$leaf_thickness)
sla_cutoff <- min(boxplot.stats(clean_leaf$sla)$out)
ldmc_cutoff <- min(boxplot.stats(clean_leaf$ldmc)$out)
leaf_thickness_cutoff <- min(boxplot.stats(clean_leaf$leaf_thickness)$out)

# filter out outliers (and scale)
range01 <- function(x){(x-min(x))/(max(x)-min(x))} # scale between 0 and 1

key_leaf_t <- clean_leaf %>%
  filter(sla < sla_cutoff) %>%
  filter(ldmc < ldmc_cutoff) %>%
  filter(leaf_thickness < leaf_thickness_cutoff) %>%
  mutate(across(all_of(c('ldmc', 'leaf_thickness', 'sla')), ~ as.numeric(range01(.x)))) # scale all vars

### find species with > 10 measurements
spp_with_10plus <- key_leaf_t %>%
  count(species) %>%
  filter(n >= 10)

key_leaf_t <- key_leaf_t %>% 
  filter(species %in% spp_with_10plus$species) # keep species with 10+ measurements

length(unique(key_leaf_t$species))

# 3. Define traits ----
traits <- c("leaf_thickness", "sla", "ldmc")

metadata <- c("id", "aspect", "site_id", "elevation_m_asl", "plant_id", "species", "family", "growth_form")

# 4. Generate ALL possible pairs----
pairs <- expand_grid(trait_x = traits, trait_y = traits) |> 
  filter(trait_x != trait_y) |>  # Remove self-pairs
  mutate(traits = paste(trait_y, trait_x, sep = "~")) |> 
  filter(traits %in% c("leaf_thickness~sla", "ldmc~sla", 'ldmc~leaf_thickness'))

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
             data = key_leaf_t,
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
nested_spp <- key_leaf_t |> 
  group_by(species) |> 
  nest()

## Expand grid to species × model grid
species_grid <- crossing(
  species = unique(key_leaf_t$species),
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

select_spp_plot_data <- spp_plot_data %>%
  filter(species %in% c("Eragrostis capensis", "Harpochloa falx", "Themeda triandra", "Helichrysum pilosellum",  "Senecio glaberrimus")) 

other_spp_plot_data <- spp_plot_data %>%
  filter(!species %in% c("Eragrostis capensis", "Harpochloa falx", "Themeda triandra", "Helichrysum pilosellum",  "Senecio glaberrimus")) 

global_plot_data <- results_best |> 
  mutate(
    data = map(predicted, ~ mutate(key_leaf_t, predicted = .x))
  ) |> 
  select(traits, significant, data) |> 
  unnest(data)

#Species order
select_spp_plot_data$species <- factor(select_spp_plot_data$species, levels = c("Eragrostis capensis", "Harpochloa falx", "Themeda triandra", "Helichrysum pilosellum",  "Senecio glaberrimus"))
select_spp_plot_data <- select_spp_plot_data %>%
  mutate(species = case_when(
    species == "Eragrostis capensis" ~ "ERCA",
    species == "Harpochloa falx" ~ "HAFA",
    species == "Themeda triandra" ~ "THTR",
    species == "Helichrysum pilosellum" ~ "HEPI",
    species == "Senecio glaberrimus" ~ "SEGL"
  ))

# Create nice labels (customize as needed)
label_lookup <- c(
  "sla" = "SLA (cm² g^-1)",
  "ldmc" = "LDMC (mg g^-1)",
  "leaf_thickness" = "Leaf thickness (mm)"
)

result_plots <- list()

for(i in 1:length(pairs$traits)){
  trait_pair <- pairs$traits[[i]]
  x_var <- str_split(trait_pair, '~')[[1]][2]
  y_var <- str_split(trait_pair, '~')[[1]][1]
  
  x_label <- label_lookup[x_var]
  y_label <- label_lookup[y_var]
  
  other_spp_data <- other_spp_plot_data |> 
    filter(traits == trait_pair)
  
  select_spp_data <- select_spp_plot_data |> 
    filter(traits == trait_pair)
  
  global_data <- global_plot_data |> 
    filter(traits == trait_pair)
  
  plot <- ggplot(other_spp_data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    ## Raw data
    geom_point(color = 'gray',alpha = 0.1,size = 2) +
    ## Species-level predictions
    geom_line(aes(y = predicted, linetype = significant, group = species), color = 'gray', alpha = 0.8, linewidth = 1) +
    ## Selected species predictions
    geom_line(data = select_spp_data, aes(y = predicted,color = species,linetype = significant),linewidth = 1) +
    ## Global prediction
    geom_line(data = global_data,aes(x = .data[[x_var]],y = predicted,linetype = significant), linewidth = 1.2, color = "black") +
    scale_color_manual(values = c("#42049EFF","#8204A7FF","#B6308BFF","#F79143FF","#FCCE25FF"),
                       name = "Species") +
    scale_linetype_manual(values = c("TRUE" = 1, "FALSE" = 2), guide = 'none') +
    labs(x = x_label, y = y_label) +
    theme_classic()
  
  result_plots[[i]] <- plot
}

combo_plot <- result_plots[[2]] + result_plots[[3]] + result_plots[[1]] + 
  plot_layout(guides = 'collect') &
  plot_annotation(tag_levels = 'A', tag_suffix = '.')
combo_plot

ggsave(filename = 'mixed_models_trait_trait_all_leaves.png',
       path = 'results/img/pairwise/',
       width = 9.59, height = 2.8, dpi = 320)

#### Tidy up model summaries ----
model_outputs <- modelsummary(results_best$model, 
                              estimate = "{estimate} [{conf.low}, {conf.high}] {stars}",
                              output = 'data.frame')
names(model_outputs)[4:6] <- results_best$model_name

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
  filter(!is.na(value)) %>%
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
  select(term, 6, 5, 4)

# Add letters on that match to the figure
names(model_outputs_clean)[2:4] <- paste0(LETTERS[1:3], '. ', names(model_outputs_clean)[2:4])

# Transpose format
model_outputs_wide <- model_outputs_clean %>%
  pivot_longer(cols = 2:4, names_to = 'model') %>%
  pivot_wider(names_from = 'term', values_from = 'value')

writexl::write_xlsx(model_outputs_wide, 'results/tab/trait_trait_mixed_ALL_LEAVES_model_outputs_raw_GLOBAL.xlsx')

#### b) Species model summaries ----
spp_model_outputs <- modelsummary(results_best_spp$model, 
                                  estimate = "{estimate} [{conf.low}, {conf.high}] {stars}",
                                  output = 'data.frame')
names(spp_model_outputs)[4:198] <- paste(results_best_spp$model_name, results_best_spp$species)

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

# create indices for selecting columns
x <- c(3+3, 3+2, 3+1)
result <- rep(x, times = 65) + rep(3 * (0:64), each = length(x))

spp_model_outputs_clean <- spp_model_outputs_clean %>%
  select(term, all_of(result))

# Add letters on that match to the figure
names(spp_model_outputs_clean)[2:196] <- paste0(rep(LETTERS[1:3], 65), '. ', names(spp_model_outputs_clean)[2:196])

# Transpose format
spp_model_outputs_wide <- spp_model_outputs_clean %>%
  pivot_longer(cols = 2:196, names_to = 'model') %>%
  pivot_wider(names_from = 'term', values_from = 'value') %>%
  mutate(species = str_extract(model, "(?<=\\)).*"),
         model = str_remove(model, "(?<=\\)).*")) %>%
  dplyr::select(1, 11, 1:10) %>%
  arrange(model)

writexl::write_xlsx(spp_model_outputs_wide, 'results/tab/trait_trait_mixed_ALL_LEAVES_model_outputs_raw_SPECIES.xlsx')
 