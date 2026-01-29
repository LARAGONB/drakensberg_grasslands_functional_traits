#Since there is measurement error in all trait measurements (Warton et al., 2006), we used standardised major axis regression (SMA) models in the lmodel2 package (Legendre, 2024). We modelled both the inter- (i.e. global patterns) and intra-specific (i.e. species-level) covariation between trait pairs using 10,000 permutations.

library(patchwork)
library(tidyverse)
library(lmodel2)
library(ggplot2)
source('scripts/functions/run_sma_multitraits.R')

#### load in the trait data
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

#### PART 1: Run for main figure ----
# select out the trait pair combinations we want to test:
# 1 leaves, 2 roots, 3 analogous leaf-root combinations
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

# run the SMA test
results <- run_sma_multitraits(
  data = trait_data_wide,
  species_col = "species",
  trait_pairs = trait_list,
  nperm = 1
)

# access the results table
results$stats_table

# save the model summaries
# write_csv(results$stats_table, 'results/tab/sma_trait_trait_relationships.csv')

# make a combined plot of all SMAs
combo_plot <- results$plots[[1]] + results$plots[[2]] + results$plots[[3]] + 
  results$plots[[4]] + results$plots[[5]] + results$plots[[6]] +
  results$plots[[7]] + results$plots[[8]] + results$plots[[9]] +
  plot_layout(guides = 'collect') &
  plot_annotation(tag_levels = 'A', tag_suffix = '.')
combo_plot

ggsave(filename = 'sma_trait_trait.png',
       path = 'results/img/',
       width = 9.59, height = 8.00, dpi = 320)

#### PART 2: Run for Plant Size traits ----
PS_trait_list <- list(
  c('BG', 'AG'),
  c('RDepth', 'VHeight'),
  c('AG', 'VHeight'),
  c('BG', 'RDepth'),
  c('BG:AG', 'SLA'),
  c('BG:AG', 'SRL'),
  c('BG', 'SLA'),
  c('BG', 'SRL'),
  c('AG', 'SLA'),
  c('AG', 'SRL')
)

# run the SMA test
PS_results <- run_sma_multitraits(
  data = filter(trait_data_wide, AG < 40),
  species_col = "species",
  trait_pairs = PS_trait_list,
  nperm = 10
)

# access the results table
PS_results$stats_table
# save the model summaries
write_csv(PS_results$stats_table, 'results/tab/sma_PLANT_SIZE_trait_trait_relationships.csv')

# make a combined plot of all SMAs
PS_legend <- cowplot::get_legend(PS_results$plots[[1]])

PS_combo_plot <- PS_results$plots[[1]] + cowplot::ggdraw(PS_legend) +
  PS_results$plots[[2]] + PS_results$plots[[3]] + 
  PS_results$plots[[6]] + PS_results$plots[[7]] +
  PS_results$plots[[4]] + PS_results$plots[[5]] + 
  # plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = list(c('A','','B', 'C', 'D', 'E', 'F', 'G')), tag_suffix = '.') +
  plot_layout(nrow = 4) &
  theme(legend.position = 'none') 
PS_combo_plot

ggsave(filename = 'sma_trait_trait_PLANT_SIZE.png',
       path = 'results/img/',
       width = 6.5, height = 10, dpi = 320)
