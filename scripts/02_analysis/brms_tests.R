### 1. Set up ----
library(tidyverse)
library(osfr)
library(lmodel2)
library(patchwork)
# install.packages('brms')
# install.packages('tidybayes')
library(brms)
library(tidybayes)

### 2. Load and format data ----
### retrieve raw data files from the OSF project page
osf_retrieve_node('hk2cy') %>%
  osf_ls_files(path = 'iv_aboveground_traits/') %>%
  filter(!str_detect(name, 'experiment')) %>% # don't download the raw scan files
  osf_download(path = 'data/raw/all_leaf_traits/', conflicts = 'overwrite')

### load in leaf trait field data
leaf <- read_csv('data/raw/all_leaf_traits/iv_PFTC7_clean_elevationgradient_traits_2023.csv')
names(leaf)
table(leaf$traits)
table(leaf$species)

### filter to select traits 
key_leaf_t <- leaf %>%
  filter(traits %in% c('ldmc', 'leaf_thickness', 'sla')) %>%
  dplyr::select(-unit) %>%
  pivot_wider(names_from = 'traits', values_from = 'value') %>%
  rename(SLA = sla, LT = leaf_thickness, LDMC = ldmc) %>% # tidy names
  filter(is.na(problem_flag)) %>%
  drop_na(c(SLA, LT, LDMC)) %>% # drop nas
  filter(LDMC < 0.8 & LDMC > 0) %>%
  filter(SLA < 1000)

### find species with > 20 measurements
spp_with_10plus <- key_leaf_t %>%
  count(species) %>%
  filter(n >= 10)

key_leaf_t <- key_leaf_t %>% 
  filter(species %in% spp_with_10plus$species) # keep species with 10+ measurements

length(unique(key_leaf_t$species))

#### load in the full trait data
trait_data_wide <- read_csv("data/processed/v_PFCT7_clean_functional_traits_2023_wide.csv") %>%
  rename(
    LT = leaf_thickness, 
    LDMC = ldmc,
    SLA = sla,
    BI = bi,
    RD = rd,
    RDMC = rdmc,
    RTD = rtd,
    SRL = srl,
    VHeight = veg_height,
    RDepth = root_depth,
    `BG:AG` = bgb_agb,
    BG = belowground_biomass,
    AG = aboveground_biomass)

# edit the levels
trait_data_wide$species <- factor(trait_data_wide$species, levels = c("Eragrostis capensis", "Harpochloa falx", "Themeda triandra", "Helichrysum pilosellum", "Senecio glaberrimus"))

#### Fit BRMs for main trait groups ----
source('scripts/functions/run_brms_multitraits.R')

# specify the trait list to run
trait_list <- list(
  c('SLA', 'LT'),
  c('SLA', 'LDMC'),
  c('LT', 'LDMC'),
  c('SRL', 'RD'),
  c('SRL', 'RTD'),
  c('RD', 'RTD'),
  c("SRL", "SLA"),
  c("RD", "LT"),
  c("RTD", "LDMC")
)

# run brms models
results <- run_trait_models(trait_data_wide, trait_list, group = "species", random_slope = TRUE)

# results[["LT_vs_SLA"]]$key_summary     # brms fixed-effect summary
# results[["LT_vs_SLA"]]$slope_table     # population & species slopes
# results[["LT_vs_SLA"]]$plot            # ggplot with dashed/solid lines

# make a combined plot of all brm plots
combo_plot <- results[[1]]$plot + results[[2]]$plot + results[[3]]$plot + 
  results[[4]]$plot + results[[5]]$plot + results[[6]]$plot +
  results[[7]]$plot + results[[8]]$plot + results[[9]]$plot +
  plot_layout(guides = 'collect') &
  plot_annotation(tag_levels = 'A', tag_suffix = '.')
combo_plot

ggsave(filename = 'brms_trait_trait.png',
       path = 'results/img/',
       width = 9.59, height = 8.00, dpi = 320)

#### Fit brms for only leaf traits on big dataset
trait_list <- list(
  c('SLA', 'LT'),
  c('SLA', 'LDMC'),
  c('LT', 'LDMC')
)

# run brms models
# use 2nd model that uses grayscale
leaf_results <- run_trait_models2(key_leaf_t, trait_list, group = "species", random_slope = TRUE)

test <- leaf_results[[1]]$slope_table

# make a combined plot of all brm plots
leaf_combo_plot <- leaf_results[[1]]$plot + leaf_results[[2]]$plot + leaf_results[[3]]$plot +
  plot_annotation(tag_levels = 'A', tag_suffix = '.')
leaf_combo_plot

ggsave(filename = 'brms_trait_trait_all_leaves.png',
       path = 'results/img/',
       width = 9.59, height = 2.8, dpi = 320)
