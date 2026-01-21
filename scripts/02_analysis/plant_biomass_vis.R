### 1. Set up ----
library(tidyverse)
library(ggpubr)
library(ggridges)
library(patchwork)
library(ggdist)
library(dunn.test)

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

species_labels <- c(
  "Eragrostis capensis" = "ERCA",
  "Harpochloa falx" = "HAFA",
  "Themeda triandra" = "THTR",
  "Helichrysum pilosellum" = "HEPI",
  "Senecio glaberrimus" = "SEGL")

####
# vis parameters
scl= 3
alp3 = 0.7
cols2<- c("#6F00A8", "#A72197",  "#DD5E66", "#EF7F4F") # 2000m: #FDCB26
lwd = 0.5

### 3. Density plots per trait by elevation (RP scans) ----
trait_data_wide %>%
  ggplot(aes(x = BG, y = species, fill = species)) +
  stat_dist_pointinterval(size = 10) +
  geom_density_ridges(alpha= alp3, scale = 2, linewidth=lwd) +
  scale_fill_manual(values = c("#42049EFF","#8204A7FF","#B6308BFF","#F79143FF","#FCCE25FF"), labels = species_labels, name = 'Species') +
  scale_x_continuous(expand = c(0,0)) +
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.key.size = unit(1, 'line')) +
  guides(fill = guide_legend(override.aes = list(size = 1))) +
  scale_y_discrete(limits = rev) +

trait_data_wide %>%
  ggplot(aes(x = AG, y = species, fill = species)) +
  stat_dist_pointinterval(size = 10) +
  geom_density_ridges(alpha= alp3, scale = 1, linewidth=lwd) +
  scale_fill_manual(values = c("#42049EFF","#8204A7FF","#B6308BFF","#F79143FF","#FCCE25FF"), labels = species_labels, name = 'Species') +
  scale_x_continuous(expand = c(0,0)) +
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.key.size = unit(1, 'line'))  +
  guides(fill = guide_legend(override.aes = list(size = 1))) +
  scale_y_discrete(limits = rev) +


trait_data_wide %>%
  ggplot(aes(x = `BG:AG`, y = species, fill = species)) +
  stat_dist_pointinterval(size = 10) +
  geom_density_ridges(alpha= alp3, scale = 2, linewidth=lwd) +
  scale_fill_manual(values = c("#42049EFF","#8204A7FF","#B6308BFF","#F79143FF","#FCCE25FF"), labels = species_labels, name = 'Species') +
  scale_x_continuous(expand = c(0,0)) +
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.key.size = unit(1, 'line')) +
  guides(fill = guide_legend(override.aes = list(size = 1))) +
  scale_y_discrete(limits = rev) +
  
  plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = 'A', tag_suffix = '.')

trait_data_wide %>%
  ggplot(aes(x = `BG:AG`, y = species, fill = species)) +
  stat_dist_pointinterval(size = 10) +
  geom_density_ridges(alpha= alp3, scale = 2, linewidth=lwd) +
  scale_fill_manual(values = c("#42049EFF","#8204A7FF","#B6308BFF","#F79143FF","#FCCE25FF"), labels = species_labels, name = 'Species') +
  scale_x_continuous(expand = c(0,0)) +
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.key.size = unit(1, 'line')) +
  guides(fill = guide_legend(override.aes = list(size = 1))) +
  scale_y_discrete(limits = rev) +
  facet_grid(~elevation_m_asl)


#### Test for differences ----  
kruskal.test(`BG:AG`~species, data = trait_data_wide)
dunn.test(
  x = trait_data_wide$`BG:AG`,
  g = trait_data_wide$species,
  method = "bh",
  kw = TRUE,
  list = TRUE
)

trait_data_wide %>%
  group_by(species) %>%
  summarise(median = median(`BG:AG`),
            min = min(`BG:AG`),
            max = max(`BG:AG`))

