################################################################################
# RDA analyses
################################################################################
#
# Lina Aragón
# linamaragonb@gmail.com
# October 14, 2025
#
# Description
################################################################################


# SET UP #######################################################################

# Load packages ----------------------------------------------------------------
#install.packages("devtools", "tidyverse", "ggplot2", "plotly", "vegan", "ggvegan", "ggrepel", "glue", "viridis", "fixest", "lmtest", "corrplot") #install if needed
#devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
pkgs <- c("devtools", "tidyverse", "ggplot2", "plotly", "vegan", "ggvegan", "ggrepel", "glue", "viridis",
          "fixest", "lmtest", "corrplot", "FactoMineR", "factoextra", "BiodiversityR",
          "cowplot", "pairwiseAdonis", "patchwork", "emmeans", "ggExtra", "multcomp",
          "multcompView", "RColorBrewer")
lapply(pkgs, library, character.only = TRUE)
remove(pkgs)

## Load data ----
### Trait data -----
trait_data_wide <- read_csv("data/processed/v_PFCT7_clean_functional_traits_2023_wide.csv")

### Climate data ----
toms_data <- read_csv("data/processed/summarised_tomst_microclimate.csv")

ndvi_summer_2023 <- read_csv("data/raw/PFTC7_sites_ndvi_S2_2023-24_summer.csv") |> 
  rename(mean_ndvi_sum = mean)

ndvi_general_2023 <- read_csv("data/raw/PFTC7_sites_ndvi_S2_2023.csv") |>  
  rename(mean_ndvi_gen = mean)

### Soil data ----
soil <- read_csv("data/raw/xii_PFTC7_clean_soil_2023.csv")

# PROCESSING ###################################################################
## Trait -----
trait_data_sum <- trait_data_wide |> 
  mutate(
    srl = log(srl),
    rtd = log(rtd),
    leaf_thickness = log(leaf_thickness),
    sla = log(sla),
    bgb_agb = log(bgb_agb)) |> 
  dplyr::select(id, species, family, growth_form, elevation_m_asl,  # Keep metadata
                leaf_thickness, ldmc, sla, 
                bi, rd, rdmc, rtd, srl, 
                veg_height, root_depth, bgb_agb) |> 
  rename(LT = leaf_thickness, LDMC = ldmc, SLA = sla, 
         BI = bi, RD = rd, RDMC = rdmc, RTD = rtd, SRL = srl, 
         VHeight = veg_height, RDepth = root_depth, `BG:AG` = bgb_agb) |> 
  group_by(species, elevation_m_asl) |>
  summarise(
    across(c(LT:'BG:AG'),
           ~ mean(.x, na.rm = TRUE)), .groups = "drop")

## Environment -----
### Climate -----
ndvi_table <- ndvi_summer_2023 |> 
  left_join(ndvi_general_2023, by = join_by(elevation_m_asl == elevation_m_asl)) |> 
  dplyr::select(elevation_m_asl, mean_ndvi_sum, mean_ndvi_gen)

climate_table <- toms_data |> 
  left_join(ndvi_table, by = join_by(elevation_m_asl == elevation_m_asl)) |> 
  pivot_wider(names_from = climate_variable,
              values_from = c(mean,sd)) |> 
  dplyr::select(-starts_with("sd"))

### Soil -----
soil_table <- soil |> 
  filter(aspect %in% "west") |> 
  summarise(mean_value = mean(value, na.rm = TRUE), 
            .by = c(elevation_m_asl, variable)) |> 
  pivot_wider(names_from = variable,
              values_from = mean_value)

### Environment ----
environment_table <- climate_table |> 
  left_join(soil_table, by = join_by(elevation_m_asl == elevation_m_asl))

## Full table ----
full_table <- trait_data_sum |>  
  left_join(environment_table, by = join_by(elevation_m_asl == elevation_m_asl))

# ANALYZE ######################################################################
## Collinearity ----------------------------------------------------------------
cor_mat <- cor(environment_table, method = "pearson")

corALL <- correlation::correlation(environment_table,
                                   include_factors = TRUE, method = "pearson")
corALL |> 
  filter(p < 0.1,
         r <= -0.7 | r >= 0.7)

## Scale and center variables --------------------------------------------------
### Trait data ----
trait_scaled <- decostand(full_table |> 
                            dplyr::select(LT:'BG:AG'), method = "standardize") |> 
  mutate(across(where(is.numeric), \(x) scale(x,
                                              center = TRUE,
                                              scale = TRUE)))

### Environment data ----
environmental_scaled <- decostand(full_table |> 
                                    dplyr::select(starts_with(c("elevation", "mean")),
                                                  cec, ph, sand, silt, tc, tn, tp), method = "standardize") |> 
  mutate(across(where(is.numeric), \(x) scale(x,
                                              center = TRUE,
                                              scale = TRUE)))


## Remove highly correlated variables ------------------------------------------
### Environment final ----
environmental_scaled_filtered <- environmental_scaled |> 
  dplyr::select(elevation_m_asl,
                mean_ndvi_sum)

## RDA analyses ----------------------------------------------------------------
env_spp <- cbind(environmental_scaled_filtered, species = full_table$species)

### With species ----
rda_m <- rda(trait_scaled ~ ., data = env_spp)
RsquareAdj(rda_m)
anova.cca(rda_m, permutations = 1000)
anova.cca(rda_m, permutations = 1000, by = "terms")
summary(rda_m)
vif.cca(rda_m)

### Constrained by spp ----
rda_partial <- rda(trait_scaled ~ elevation_m_asl + mean_ndvi_sum + Condition(species), data = env_spp)
RsquareAdj(rda_partial)$adj.r.squared*100
set.seed(4321)
anova.cca(rda_partial, permutations = 1000)
set.seed(4321)
anova.cca(rda_partial, permutations = 1000, by = "terms")
summary(rda_partial)
vif.cca(rda_partial)

### Soil----
rda_soil <- rda(trait_scaled ~ cec + ph + sand + silt + tc + tn + tp, 
                data = environmental_scaled)
RsquareAdj(rda_soil)$adj.r.squared*100
set.seed(4321)
anova.cca(rda_soil, permutations = 1000)
set.seed(4321)
anova.cca(rda_soil, permutations = 1000, by = "terms")
set.seed(4321)
anova.cca(rda_soil, permutations = 1000, by = "axis")
summary(rda_soil)
vif.cca(rda_soil)


#### Percentage of information explained by partial model --------------------------
percent_rda <- round(100*(summary(rda_partial)$cont$importance[2, 1:2]), 2)
eigenvals_rda <- eigenvals(rda_partial)[1:2]  # Constrained axes only
prop_explained <- eigenvals_rda / sum(eigenvals_rda) * 100


### Extract scores ------------------------------------------------------------
#### Site -----
site_scores <- scores(rda_partial, display = "sites", choices = 1:2)
# Add species labels to site scores
site_df <- as.data.frame(site_scores) |> 
  mutate(species = full_table$species,
         elevation = full_table$elevation_m_asl,
         species = fct_relevel(species, "Themeda triandra", after = 2))

#### Traits -----
trait_scores <- scores(rda_partial, display = "species", choices = 1:2)
# Trait loadings
trait_df <- as.data.frame(trait_scores) |> 
  rownames_to_column(var = "trait")

trait_df |> 
  filter(RDA1 <= -0.3 | RDA1 >= 0.3)

#### Environment ----
env_scores <- scores(rda_partial, display = "bp", choices = 1:2)
# Environmental vectors
env_df <- as.data.frame(env_scores) |> 
  rownames_to_column(var = "variable") |> 
  mutate(variable = case_when(
    variable == "elevation_m_asl" ~ "Elevation (m asl)",
    variable == "mean_ndvi_sum" ~ "Summer NDVI"
  ))

#### Correlate traits with RDA axes -------
trait_env_cor <- cor(trait_scaled, site_scores)

cat("Correlation between traits and RDA axes:\n")
print(round(trait_env_cor, 3))

# Which traits are most strongly correlated with environmental gradients?
trait_cor_df <- as.data.frame(trait_env_cor)
trait_cor_df$trait <- rownames(trait_cor_df)
trait_cor_df$max_cor <- apply(abs(trait_cor_df[, 1:2]), 1, max)

# VISUALIZE ####################################################################
cor_plot <- ggplot(as_tibble(cor_mat, rownames = "Var1") |> 
                     pivot_longer(!Var1,
                                  names_to = "Var2", values_to = "value"), 
                   aes(Var2, Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(colors = rev(brewer.pal(11, "PiYG")),
                       values = scales::rescale(c(-1, 0, 1)), # ensures midpoint is 0
                       limits = c(-1, 1), na.value = "grey90",
                       guide = guide_colorbar(barwidth = unit(0.4, "cm"),
                                              barheight = unit(7, "cm"))) +
  coord_equal(expand = FALSE) +
  labs(x = "",
       y = "") +
  theme_minimal() +
  theme(
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.25, size = 10,
                               color = "black"))

cor_plot

rda_plot <- ggplot(data = site_df, aes(x = RDA1, y = RDA2)) +
  geom_hline(yintercept = 0, color = "grey90", linewidth = 1) +     
  geom_vline(xintercept = 0, color = "grey90", linewidth = 1) +
  stat_ellipse(aes(group = factor(elevation), linetype = factor(elevation)),
               color = "gray50", 
               linewidth = 0.8, 
               level = 0.95) +
  geom_point(aes(shape = as.factor(elevation), colour = species), 
             size = 4) +
  scale_colour_manual(values = species_colors, name = "Species",
                      labels = species_labels) +
  scale_shape_manual(values = elevation_shapes, name = "Elevation (m asl)",
                     guide = guide_legend(
                       override.aes = list(color = "gray50"))) +
  scale_linetype_manual(values = elevation_ellipses, name = "Elevation (m asl)") +
  geom_segment(data = trait_df,
               aes(x = 0, y = 0, xend = RDA1*2.5, yend = RDA2*2.5),
               arrow = arrow(length = unit(0.15, "cm")),
               color = "grey10", alpha = 0.6) +
  geom_text_repel(data = trait_df,
                  aes(x = RDA1*3, y = RDA2*3, label = trait),
                  size = 3, color = "grey10",
                  max.overlaps = 20) +
  geom_segment(data = env_df,
               aes(x = 0, y = 0, xend = RDA1*3.5, yend = RDA2*3.5),
               arrow = arrow(length = unit(0.3, "cm")),
               color = "black", linewidth = 1) +
  geom_text_repel(data = env_df,
                  aes(x = RDA1*4.1, y = RDA2*4.1, label = variable),
                  size = 4, fontface = "bold", color = "black") +
  labs(
    x = glue("PCA1 ({round(prop_explained[1], 1)}% of constrained variance)"),
    y = glue("PCA2 ({round(prop_explained[2], 1)}% of constrained variance)")) +
  theme_bw(base_size = 16) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.ticks = element_line(linewidth = 1),
    axis.line = element_line(linewidth = 1, colour = "black"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 16),
    legend.key.width = unit(2, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(0.3, 0.3, 0.3, 0.3, "mm"),
    aspect.ratio = 1)
rda_plot

## EXPORT #######################################################################
ggsave("results/img/rda/cor_plot.png", cor_plot, 
       width = 8, height = 8, dpi = 300)

ggsave("results/img/rda/rda_plot.png", rda_plot, 
       width = 8, height = 8, dpi = 300)

ggplot 