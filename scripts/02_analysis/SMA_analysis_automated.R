library(patchwork)
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
    `BG:AG` = bgb_agb)

# edit the levels
trait_data_wide$species <- factor(trait_data_wide$species, levels = c("Eragrostis capensis", "Harpochloa falx", "Themeda triandra", "Helichrysum pilosellum", "Senecio glaberrimus"))

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
  c("RD", "SLA"),
  c("RTD", "LDMC")
)

# run the SMA test
results <- run_sma_multitraits(
  data = trait_data_wide,
  species_col = "species",
  trait_pairs = trait_list,
  nperm = 4999
)

# access the results table
results$stats_table
# save the model summaries
write_csv(results$stats_table, 'results/tab/sma_trait_trait_relationships.csv')

# make a combined plot of all SMAs
combo_plot <- results$plots[[1]] + results$plots[[2]] + results$plots[[3]] + 
  results$plots[[4]] + results$plots[[5]] + results$plots[[6]] +
  results$plots[[7]] + results$plots[[8]] + results$plots[[9]] +
  plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = 'A', tag_suffix = '.')
combo_plot

ggsave(filename = 'sma_trait_trait.png',
       path = 'results/img/',
       width = 9.59, height = 8.00, dpi = 320)
