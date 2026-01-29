### 1. Set up ----
library(tidyverse)
library(osfr)
library(lmodel2)
library(patchwork)
# install.packages('brms')
# install.packages('tidybayes')
library(brms)
library(tidybayes)
library(report)
library(sjPlot)
library(gtsummary)

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

#### Part 1: Fit BRMs for main trait groups ----
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

trait_list <- list(
  c('l_SRL', 'l_RD')
)
trait_data_wide <- trait_data_wide %>%
  mutate(l_SRL = log(SRL),
         l_RD = log(RD))

# run brms models
results <- run_trait_models(trait_data_wide, trait_list, group = "species", random_slope = TRUE)

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

trait_trait_slope_tables <- bind_rows(results[[1]]$slope_table, results[[2]]$slope_table, results[[3]]$slope_table, results[[4]]$slope_table, results[[5]]$slope_table, results[[6]]$slope_table, results[[7]]$slope_table, results[[8]]$slope_table, results[[9]]$slope_table)
trait_trait_key_summary <- bind_rows(results[[1]]$key_summary, results[[2]]$key_summary, results[[3]]$key_summary, results[[4]]$key_summary, results[[5]]$key_summary, results[[6]]$key_summary, results[[7]]$key_summary, results[[8]]$key_summary, results[[9]]$key_summary)

write_csv(trait_trait_slope_tables, 'results/tab/brms_trait_trait_slope_tables.csv')
write_csv(trait_trait_key_summary, 'results/tab/brms_trait_trait_key_summary.csv')

#### report
r <- report(results[[1]]$fit)
r
summary(results[[1]]$fit)
summary(r)
report_random(r)
report_text(r)

tab_model(results[[1]]$fit)
tbl_regression(results[[1]]$fit, tidy_fun = broom.mixed::tidy)

gtsummary::tbl_summary(results[[1]]$fit)


#### Part 2: Fit brms for only leaf traits on big dataset ----
trait_list <- list(
  c('SLA', 'LT'),
  c('SLA', 'LDMC'),
  c('LT', 'LDMC')
)

# run brms models
# use 2nd model that uses grayscale
leaf_results <- run_trait_models2(key_leaf_t, trait_list, group = "species", random_slope = TRUE)

# make a combined plot of all brm plots
leaf_combo_plot <- leaf_results[[1]]$plot + leaf_results[[2]]$plot + leaf_results[[3]]$plot +
  plot_annotation(tag_levels = 'A', tag_suffix = '.')
leaf_combo_plot

ggsave(filename = 'brms_trait_trait_all_leaves.png',
       path = 'results/img/',
       width = 9.59, height = 2.8, dpi = 320)

#### save leaf tables
leaf_slope_tables <- bind_rows(leaf_results[[1]]$slope_table, leaf_results[[2]]$slope_table, leaf_results[[3]]$slope_table)
leaf_key_summary <- bind_rows(leaf_results[[1]]$key_summary, leaf_results[[2]]$key_summary, leaf_results[[3]]$key_summary)

write_csv(leaf_slope_tables, 'results/tab/brms_leaf_slope_tables.csv')
write_csv(leaf_key_summary, 'results/tab/brms_leaf_key_summary.csv')

#### PART 3: Run for Plant Size traits ----
PS_trait_list <- list(
  c('BG', 'AG'),
  c('RDepth', 'VHeight'),
  c('AG', 'VHeight'),
  c('BG', 'RDepth'),
  c('BG_AG', 'SLA'),
  c('BG_AG', 'SRL'),
  c('BG', 'SLA'),
  c('BG', 'SRL'),
  c('AG', 'SLA'),
  c('AG', 'SRL')
)

# PS_trait_list_1 <- list(
#   c('BG', 'AG'),
#   c('RDepth', 'VHeight'),
#   c('AG', 'VHeight')
#   )
# 
# PS_trait_list_2 <- list(
#   c('BG', 'RDepth'),
#   c('BG_AG', 'SLA'),
#   c('BG_AG', 'SRL')
# )
# 
# PS_trait_list_3 <- list(
#   c('BG', 'SLA'),
#   c('BG', 'SRL'),
#   c('AG', 'SLA'),
#   c('AG', 'SRL')
# )

# filter out large AG measure
trait_data_wide_filt <- filter(trait_data_wide, AG < 40) %>%
  rename(BG_AG = `BG:AG`)

# run brms models
ps_results <- run_trait_models(trait_data_wide_filt, PS_trait_list, group = "species", random_slope = TRUE)

# make a combined plot of all brm plots
PS_legend <- cowplot::get_legend(ps_results[[1]]$plot)

ps_combo_plot <- ps_results[[1]]$plot + ps_results[[2]]$plot + cowplot::ggdraw(PS_legend) + ps_results[[4]]$plot + ps_results[[8]]$plot + ps_results[[7]]$plot + 
  ps_results[[3]]$plot + ps_results[[10]]$plot + ps_results[[9]]$plot +
  # plot_layout(guides = 'collect') &
  plot_annotation(tag_levels = list(c('A.','B.','', 'C.', 'D.', 'E.', 'F.', 'G.','H.'))) &
  theme(legend.position = 'none') 
ps_combo_plot

ggsave(filename = 'brms_ps_trait_trait.png',
       path = 'results/img/',
       width = 9.59, height = 8.00, dpi = 320)

#### save ps tables
ps_slope_tables <- bind_rows(ps_results[[1]]$slope_table , ps_results[[2]]$slope_table , ps_results[[4]]$slope_table , ps_results[[8]]$slope_table , ps_results[[7]]$slope_table , 
                               ps_results[[3]]$slope_table , ps_results[[10]]$slope_table , ps_results[[9]]$slope_table)
ps_key_summary <- bind_rows(ps_results[[1]]$key_summary , ps_results[[2]]$key_summary ,  ps_results[[4]]$key_summary , ps_results[[8]]$key_summary , ps_results[[7]]$key_summary , 
                              ps_results[[3]]$key_summary , ps_results[[10]]$key_summary , ps_results[[9]]$key_summary)

write_csv(ps_slope_tables, 'results/tab/brms_ps_slope_tables.csv')
write_csv(ps_key_summary, 'results/tab/brms_ps_key_summary.csv')

