### 1. Set up ----
library(tidyverse)
library(osfr)
library(lmodel2)
library(patchwork)
source('scripts/functions/run_sma_multitraits.R')

### retrieve raw data files from the OSF project page
osf_retrieve_node('hk2cy') %>%
  osf_ls_files(path = 'iv_aboveground_traits/') %>%
  filter(!str_detect(name, 'experiment')) %>% # don't download the raw scan files
  osf_download(path = 'data/raw/all_leaf_traits/', conflicts = 'overwrite')

### 2. Format data ----
### load in leaf trait field data
leaf <- read_csv('data/raw/all_leaf_traits/iv_PFTC7_clean_elevationgradient_traits_2023.csv')
names(leaf)
table(leaf$traits)
table(leaf$species)

### find species with > 10 measurements
spp_with_10plus <- leaf %>%
  count(species) %>%
  filter(n >= 20)

### filter to select traits 
key_leaf_t <- leaf %>%
  filter(traits %in% c('ldmc', 'leaf_thickness', 'sla')) %>%
  dplyr::select(-unit) %>%
  pivot_wider(names_from = 'traits', values_from = 'value') %>%
  rename(SLA = sla, LT = leaf_thickness, LDMC = ldmc) %>% # tidy names
  filter(is.na(problem_flag)) %>%
  drop_na(c(SLA, LT, LDMC)) %>% # drop nas
  filter(species %in% spp_with_10plus$species) %>% # keep species with 10+ measurements
  filter(LDMC < 0.8 & LDMC > 0) %>%
  filter(SLA < 1000)

length(unique(key_leaf_t$species))
#### 3. Run for SMA ----
# select out the trait pair combinations we want to test:
trait_list <- list(
  c('SLA', 'LT'),
  c('SLA', 'LDMC'),
  c('LT', 'LDMC')
)

# run the SMA test
source('scripts/functions/run_sma_multitraits.R')
results <- run_sma_multitraits2(
  data = key_leaf_t,
  species_col = "species",
  trait_pairs = trait_list,
  nperm = 1
)

# access the results table
# results$stats_table

#
results$plots[[1]] +
results$plots[[2]] +
results$plots[[3]] +
  plot_annotation(tag_levels = 'A', tag_suffix = '.')

ggsave(filename = 'sma_trait_trait_all_leaves.png',
       path = 'results/img/',
       width = 9.59, height = 2.80, dpi = 320)

# #### test this dataset against Vile (2005) predictions
# key_leaf_t_vile <- key_leaf_t %>%
#   mutate(inv_prod_sla_x_ldmc = 1/(SLA*LDMC),
#          log_LT = log(1+LT),
#          log_inv_prod = log(1+inv_prod_sla_x_ldmc))
# 
# key_leaf_t_vile %>%
#   filter(inv_prod_sla_x_ldmc < 0.15) %>%
#   ggplot(aes(x = inv_prod_sla_x_ldmc, y = LT)) +
#   geom_point()
# 
# key_leaf_t_vile %>%
#   filter(inv_prod_sla_x_ldmc < 0.15) %>%
#   ggplot(aes(x = log_inv_prod, y = log_LT)) +
#   geom_point()

