################################################################################
# PCA analyses
################################################################################
#
# Lina Aragón
# linamaragonb@gmail.com
# September 24, 2025
#
# Description
################################################################################

# 1. Load libraries ----

#install.packages("devtools", "tidyverse", "ggplot2", "plotly", "vegan", "ggvegan", "ggrepel", "glue", "viridis", "fixest", "lmtest", "corrplot") #install if needed
# devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
pkgs <- c("devtools", "tidyverse", "ggplot2", "plotly", "vegan", "ggvegan", "ggrepel", "glue", "viridis",
          "fixest", "lmtest", "corrplot", "FactoMineR", "factoextra", "BiodiversityR",
          "cowplot", "pairwiseAdonis")
lapply(pkgs, library, character.only = TRUE)
remove(pkgs)

# 2. Load data ----
trait_data <- read_csv("data/processed/v_PFCT7_clean_functional_traits_2023.csv")
trait_data_wide <- read_csv("data/processed/v_PFCT7_clean_functional_traits_2023_wide.csv")

# 3. PCA all traits ----
## PCA using vegan (rda with no constraints = PCA) ----
## Remember to scale to unit variance
pca_output <- trait_data_wide |> 
  select(root_depth, veg_height, rd, bi, srl, rtd, rdmc, sla, ldmc, leaf_thickness, bgb_agb) |> 
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
    `BG:AG` = bgb_agb) |> 
  rda(scale = TRUE)
summary(pca_output)

## How many PCs should we keep? ----
### Extract eigenvalues and their percentage
ev <- eigenvals(pca_output)
e_B <- eigenvals(pca_output)/sum(eigenvals(pca_output))
### Extract PCs that explain 90% of the variation in the data
k <- which(cumsum(ev) / sum(ev) >= 0.9)[1]

## Tables for ind. loadings and traits loadings ----
## Create a table containing the loading of each individual (sites) in each PC

pca_sites <- as_tibble(bind_cols(
  trait_data_wide |> 
    mutate(species = fct_relevel(species, "Themeda triandra", after = 2)) |> 
    select(id, species, family, growth_form, elevation_m_asl),
  scores(pca_output, display = "sites", choices = 1:k, scaling = 2)))

pca_traits <- scores(pca_output, display = "species", choices = 1:k, scaling = 2) |> 
  as.data.frame() |> 
  rownames_to_column(var = "traits") 

## PCA significance
PCAsignificance(pca_output)
plot1 <- ordiplot(pca_output, choices = c(1,2), scaling = 1)
ordiequilibriumcircle(pca_output, plot1) #Which traits are more important in each PC


## PERMANOVA  ----
### Matrix of distances ----
all_distance <- dist(pca_sites |> 
                       select(starts_with("PC")))
### Species alone ----
set.seed(1)
adonis2(all_distance ~ species, data = pca_sites, permutations = 4999)

# Multivarite homogeneity
all_sp_bd <- betadisper(all_distance, pca_sites$species)
anova(all_sp_bd)
permutest(all_sp_bd, 999)
plot(all_sp_bd)

# Plot showing differences among spp
score_spp <- scores(pca_output, display = "sites", scaling = 1, choices = 1:2)
plot(score_spp, type = "n") 
points(score_spp, col = as.integer(pca_sites$species), pch = 19)
ordiellipse(score_spp, pca_sites$species, kind = "se", conf = 0.95, draw = "polygon",
            col = 1:5, border = 1:5, label = TRUE)

#### Pairwise comparisons among species ----
set.seed(1)
pairwise.adonis2(all_distance ~ species, 
                 data = as.data.frame(pca_sites),
                 permutations = 4999, 
                 p.adjust.m = "holm")

### Elevation alone ----
set.seed(1)
adonis2(all_distance ~ elevation_m_asl, data = pca_sites, permutations = 4999)

# Multivarite homogeneity
all_ele_bd <- betadisper(all_distance, pca_sites$elevation_m_asl)
perm_all_ele_bd  <- permutest(all_ele_bd, permutations = 4999)
anova(all_ele_bd)
permutest(all_ele_bd, 999)
plot(all_ele_bd)

# Plot showing differences among elevations
score_ele <- scores(pca_output, display = "sites", scaling = 1, choices = 1:2)
plot(score_ele, type = "n") 
points(score_ele, col = as.integer(pca_sites$elevation_m_asl), pch = 19)
ordiellipse(score_ele, pca_sites$elevation_m_asl, kind = "se", conf = 0.95, draw = "polygon",
            col = 1:5, border = 1:5, label = TRUE)

#### Pairwise comparisons among elevations ----
set.seed(1)
pairwise.adonis2(all_distance ~ elevation_m_asl, 
                 data = as.data.frame(pca_sites),
                 permutations = 4999, 
                 p.adjust.m = "holm")

### Marginal differences species + elevation ----
set.seed(1)
adonis2(all_distance ~ species + elevation_m_asl, data = as.data.frame(pca_sites), 
        permutations = 4999,
        by = "margin")


pairwise.adonis2(all_distance ~ species + elevation_m_asl, 
                 data = as.data.frame(pca_sites),
                 p.adjust.m = "holm")


### Species pooled “within-elevation” test ----
set.seed(1)
adonis2(all_distance ~ species, data = as.data.frame(pca_sites), permutations = 4999,
        strata = pca_sites$elevation_m_asl)

#### Pairwise comparisons among species ----
set.seed(1)
pairwise.adonis2(all_distance ~ species, 
                 data = as.data.frame(pca_sites),
                 strata = "elevation_m_asl", 
                 p.adjust.m = "holm")


### Elevation pooled “within-species” test ----
set.seed(1)
adonis2(all_distance ~ elevation_m_asl, data = as.data.frame(pca_sites), 
        permutations = 4999,
        strata = pca_sites$species)

#### Pairwise comparisons among elevations within species ----
set.seed(1)
pairwise.adonis2(all_distance ~ elevation_m_asl, 
                 data = as.data.frame(pca_sites),
                 strata = "species", 
                 p.adjust.m = "holm")

## PCA Plot ----
xlim_equal <- c(pca_sites$PC2, pca_sites$PC3) |> 
  abs() |> 
  max(na.rm = TRUE) |> 
  (\(m) c(-m, m))()

### PC1 & PC2 ----
pca12_all <- pca_sites |> 
  ggplot(aes(x = PC1, y = PC2, 
             colour = species)) +
  geom_point(aes(shape = factor(elevation_m_asl)), size = 4) +
  scale_shape_manual(values = elevation_shapes, name = "Elevation (m asl)") +  # Adjust the number of shapes to match your elevation count
  stat_ellipse(aes(group = species, colour = species), size = 0.8) +
  geom_segment(data = pca_traits,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.5, "cm")),
               size = 1,
               colour = "grey20",
               inherit.aes = FALSE) +
  geom_text_repel(data = pca_traits,
                  aes(x = PC1 * 1.1, y = PC2 * 1.1, label = traits),
                  size = 4,
                  fontface = "bold",
                  inherit.aes = FALSE, 
                  colour = "black") +
  coord_equal() +
  scale_colour_manual(values = species_colors, name = "Species",
                      labels = species_labels) +
  labs(x = glue("PCA1 ({round(e_B[1] * 100, 1)}%)"),
       y = glue("PCA2 ({round(e_B[2] * 100, 1)}%)")) +
  theme_bw(base_size = 16) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.ticks = element_line(linewidth = 1),
    axis.line = element_line(linewidth = 1, colour = "black"),
    legend.text = element_text(size = 16),
    plot.margin = margin(2,2,2,2),
    aspect.ratio = 1) +
  coord_cartesian(xlim = xlim_equal)
 
pca12_all

#### Save plots ----
# ggsave("results/img/pca12_all_tiff.tiff", pca12_all,
#        width = 20, height = 20, units = "cm", dpi = 300)
ggsave("results/img/pca12_all_png.png", pca12_all,
       width = 20, height = 20, units = "cm", dpi = 300)

### PC1 & PC3 ----
pca13_all <- pca_sites |> 
  ggplot(aes(x = PC1, y = PC3, 
             colour = species)) +
  geom_point(aes(shape = factor(elevation_m_asl)), size = 4) +
  scale_shape_manual(values = elevation_shapes, name = "Elevation (m asl)") +  # Adjust the number of shapes to match your elevation count
  stat_ellipse(aes(group = species, colour = species), size = 0.8) +
  geom_segment(data = pca_traits,
               aes(x = 0, y = 0, xend = PC1, yend = PC3),
               arrow = arrow(length = unit(0.5, "cm")),
               size = 1,
               colour = "grey20",
               inherit.aes = FALSE) +
  geom_text_repel(data = pca_traits,
                  aes(x = PC1 * 1.1, y = PC3 * 1.1, label = traits),
                  size = 4,
                  fontface = "bold",
                  inherit.aes = FALSE, 
                  colour = "black") +
  coord_equal() +
  scale_colour_manual(values = species_colors, name = "Species",
                      labels = species_labels) +
  labs(x = glue("PCA1 ({round(e_B[1] * 100, 1)}%)"),
       y = glue("PCA3 ({round(e_B[3] * 100, 1)}%)")) +
  theme_bw(base_size = 16) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.ticks = element_line(linewidth = 1),
    axis.line = element_line(linewidth = 1, colour = "black"),
    legend.text = element_text(size = 16),
    plot.margin = margin(2,2,2,2),
    aspect.ratio = 1) +
  coord_cartesian(xlim = xlim_equal)

pca13_all

#### Save plots ----
# ggsave("results/img/pca13_all_tiff.tiff", pca13_all,
#        width = 20, height = 20, units = "cm", dpi = 300)
ggsave("results/img/pca13_all_png.png", pca13_all,
       width = 20, height = 20, units = "cm", dpi = 300)

### PC2 & PC3 ----
pca23_all <- pca_sites |> 
  ggplot(aes(x = PC2, y = PC3, 
             colour = species)) +
  geom_point(aes(shape = factor(elevation_m_asl)), size = 4) +
  scale_shape_manual(values = elevation_shapes, name = "Elevation (m asl)") +  # Adjust the number of shapes to match your elevation count
  stat_ellipse(aes(group = species, colour = species), size = 0.8) +
  geom_segment(data = pca_traits,
               aes(x = 0, y = 0, xend = PC2, yend = PC3),
               arrow = arrow(length = unit(0.5, "cm")),
               size = 1,
               colour = "grey20",
               inherit.aes = FALSE) +
  geom_text_repel(data = pca_traits,
                  aes(x = PC2 * 1.1, y = PC3 * 1.1, label = traits),
                  size = 4,
                  fontface = "bold",
                  inherit.aes = FALSE, 
                  colour = "black") +
  coord_equal() +
  scale_colour_manual(values = species_colors, name = "Species",
                      labels = species_labels) +
  labs(x = glue("PCA2 ({round(e_B[2] * 100, 1)}%)"),
       y = glue("PCA3 ({round(e_B[3] * 100, 1)}%)")) +
  theme_bw(base_size = 16) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.ticks = element_line(linewidth = 1),
    axis.line = element_line(linewidth = 1, colour = "black"),
    legend.text = element_text(size = 16),
    plot.margin = margin(2,2,2,2),
    aspect.ratio = 1) +
  coord_cartesian(xlim = xlim_equal)

pca23_all

#### Save plots ----
# ggsave("results/img/pca23_all_tiff.tiff", pca23_all,
#        width = 20, height = 20, units = "cm", dpi = 300)
ggsave("results/img/pca23_all_png.png", pca23_all,
       width = 20, height = 20, units = "cm", dpi = 300)


plots_pca_all <- c(pca12_all, pca13_all, pca23_all)
pca_all_legend <- get_legend(plots_pca_all[[1]] + theme(legend.position = "right"))
plots_pca_all_nolegend <- lapply(plots_pca_all, function(p) p + theme(legend.position = "none"))


plot_all_pca <- plot_grid(plots_pca_all_nolegend[[1]], plots_pca_all_nolegend[[2]],
                          plots_pca_all_nolegend[[3]], pca_all_legend,
                          ncol = 2,
                          rel_widths = c(1, 1), 
                          rel_heights = c(1, 1),
                          labels = c("A.", "B.", "C.", ""),
                          align = "hv",
                          axis = "tblr") + 
  theme(plot.background = element_rect(fill = "white", colour = NA))

plot_all_pca 

#### Save all plots ----
ggsave("results/img/plot_all_pca_tiff.tiff", plot_all_pca,
       width = 25, height = 25, units = "cm", dpi = 300)
ggsave("results/img/plot_all_pca_png.png", plot_all_pca,
       width = 25, height = 25, units = "cm", dpi = 300)

# 4. PCA Roots ----
## PCA using vegan (rda with no constraints = PCA) ----
## Remember to scale to unit variance
roots_traits <- trait_data_wide |> 
  select(rd, bi, srl, rtd, rdmc) |> 
  rename(BI = bi, RD = rd, RDMC = rdmc, RTD = rtd, SRL = srl)

roots_sum <- roots_traits |> 
  pivot_longer(everything(), names_to = "trait", values_to = "value") |> 
  group_by(trait) |> 
  summarise(min = min(value),
            max = max(value),
            mean = mean(value),
            median = median(value),
            zero = sum(value == 0),
            skewness = e1071::skewness(value)) ; roots_sum

roots_trans <- roots_traits %>%
  mutate(
    SRL_log = log(SRL),
    RTD_log = log(RTD)) |> 
  select(-SRL, -RTD)

roots_sum_trans <- roots_trans |> 
  pivot_longer(everything(), names_to = "trait", values_to = "value") |> 
  group_by(trait) |> 
  summarise(min = min(value),
            max = max(value),
            mean = mean(value),
            median = median(value),
            zero = sum(value == 0),
            skewness = e1071::skewness(value)); roots_sum_trans

pca_roots <- roots_trans |> 
  rda(scale = TRUE, center = TRUE)

summary(pca_roots)

## How many PCs should we keep? ----
### Extract eigenvalues and their percentage
ev_roots <- eigenvals(pca_roots); ev_roots
e_B_roots <- eigenvals(pca_roots)/sum(eigenvals(pca_roots)); e_B_roots
### Extract PCs that explain 90% of the variation in the data
k_roots <- which(cumsum(ev_roots) / sum(ev_roots) >= 0.9)[1]; k_roots

## Tables for ind. loadings and traits loadings ----
## Create a table containing the loading of each invidual (sites) in each PC

pca_sites_roots <- as_tibble(bind_cols(
  trait_data_wide |> 
    mutate(species = fct_relevel(species, "Themeda triandra", after = 2)) |> 
    select(id, species, family, growth_form, elevation_m_asl),
  scores(pca_roots, display = "sites", choices = 1:k_roots, scaling = 2)))

pca_traits_roots <- scores(pca_roots, display = "species", choices = 1:k_roots, scaling = 2) |> 
  as.data.frame() |> 
  rownames_to_column(var = "traits") 

## PCA significance
sig_roots_pca <- PCAsignificance(pca_roots)
plot1_roots <- ordiplot(pca_roots, choices = c(1,2), scaling = 1)
ordiequilibriumcircle(pca_roots,plot1_roots) #Which traits are more important in each PC


barplot (sig_roots_pca[c('percentage of variance', 'broken-stick percentage'), ], beside = T, 
         xlab = 'PCA axis', ylab = 'explained variation [%]', col = c('grey', 'black'), 
         legend = TRUE)

## PERMANOVA  ----
### Matrix of distances ----
roots_distance <- dist(pca_sites_roots |> 
                       select(starts_with("PC")), method = "euclidean")

#### Multivarite homogeneity ----
# elevation
roots_ele_bd <- betadisper(roots_distance, pca_sites_roots$elevation_m_asl)
perm_roots_ele_bd  <- permutest(roots_ele_bd, permutations = 4999)
print(perm_roots_ele_bd)
plot(roots_ele_bd)

# species
roots_spp_bd <- betadisper(roots_distance, pca_sites_roots$species)
perm_roots_spp_bd  <- permutest(roots_spp_bd, permutations = 4999)
print(perm_roots_spp_bd)
plot(roots_spp_bd)

#### Permanova Marginal differences species + elevation ----
set.seed(4321)
permanova_roots <- adonis2(roots_distance ~ species + elevation_m_asl, 
        data = as.data.frame(pca_sites_roots), 
        permutations = 10000,
        by = "margin")
print(permanova_roots)


set.seed(4321)
pairwise.adonis2(roots_distance ~ species, 
                 data = as.data.frame(pca_sites_roots),
                 strata = "elevation_m_asl",
                 nperm = 10000,
                 p.adjust.m = "BH")

## PCA Plot ----
xlim_equal_roots <- c(pca_traits_roots$PC1, pca_traits_roots$PC2) |> 
  abs() |> 
  max(na.rm = TRUE) |> 
  (\(m) c(-m, m))()

### PC1 & PC2 ----
pca12_roots <- pca_sites_roots |> 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_hline(yintercept = 0, color = "grey90", linewidth = 0.5, linetype = "dashed") +     
  geom_vline(xintercept = 0, color = "grey90", linewidth = 0.5, linetype = "dashed") + 
  geom_point(aes(shape = as.factor(elevation_m_asl)), 
             size = 3, color = "grey70") +
  stat_ellipse(aes(linetype = species), 
               linewidth = 0.8, show.legend = TRUE, color = "grey30") +
  geom_segment(data = pca_traits_roots |> 
                 filter(!traits %in% "BI"),
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.18, "cm")),
               linewidth = 1.2,
               inherit.aes = FALSE,
               color = group_colors["Roots"]) +
  geom_text_repel(data = pca_traits_roots |> 
                    filter(!traits %in% "BI"),
                  aes(x = PC1, y = PC2, label = traits),
                  nudge_x = 0.4,
                  size = 5,
                  show.legend = FALSE,
                  fontface = "bold", 
                  color = group_colors["Roots"]) +
  coord_equal() +
  scale_shape_manual(values = elevation_shapes, name = "Elevation",
                    labels = elevation_labels) +
  scale_linetype_manual(values = species_ellipses, 
                        name = "Species",
                        labels = species_labels) +
  labs(x = glue("PCA1 ({round(e_B_roots[1] * 100, 1)}%)"),
       y = glue("PCA2 ({round(e_B_roots[2] * 100, 1)}%)")) +
  guides(
    shape = guide_legend(override.aes = list(linetype = "blank", 
                                             colour = "grey70")),
    linetype = guide_legend(override.aes = list(shape = NA))) +
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(linewidth = 0.5),
    axis.line = element_line(linewidth = 1, colour = "black"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(),
    legend.key.width = unit(1.2, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(1,1,1,1, "mm"),
    aspect.ratio = 1)
 

pca12_roots

#### Save plots ----
# ggsave("results/img/pca12_roots_tiff.tiff", pca12_roots,
#        width = 20, height = 20, units = "cm", dpi = 300)
ggsave("results/img/pca12_roots_png.png", pca12_roots,
       width = 15, height = 15, units = "cm", dpi = 300)

### PC1 & PC3 ----
pca13_roots <- pca_sites_roots |> 
  ggplot(aes(x = PC1, y = PC3, 
             colour = species)) +
  geom_point(aes(shape = factor(elevation_m_asl)), size = 4) +
  scale_shape_manual(values = elevation_shapes, name = "Elevation (m asl)") +  # Adjust the number of shapes to match your elevation count
  stat_ellipse(aes(group = species, colour = species), size = 0.8) +
  geom_segment(data = pca_traits_roots,
               aes(x = 0, y = 0, xend = PC1, yend = PC3),
               arrow = arrow(length = unit(0.5, "cm")),
               size = 1,
               colour = "grey20",
               inherit.aes = FALSE) +
  geom_text_repel(data = pca_traits_roots,
                  aes(x = PC1 * 1.1, y = PC3 * 1.1, label = traits),
                  size = 4,
                  fontface = "bold",
                  inherit.aes = FALSE, 
                  colour = "black") +
  coord_equal() +
  scale_colour_manual(values = species_colors, name = "Species",
                      labels = species_labels) +
  labs(x = glue("PCA1 ({round(e_B_roots[1] * 100, 1)}%)"),
       y = glue("PCA3 ({round(e_B_roots[3] * 100, 1)}%)")) +
  theme_bw(base_size = 16) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.ticks = element_line(linewidth = 1),
    axis.line = element_line(linewidth = 1, colour = "black"),
    legend.text = element_text(size = 16),
    plot.margin = margin(2,2,2,2),
    aspect.ratio = 1) +
  coord_cartesian(xlim = xlim_equal_roots)

pca13_roots

#### Save plots ----
# ggsave("results/img/pca13_roots_tiff.tiff", pca13_roots,
#        width = 20, height = 20, units = "cm", dpi = 300)
ggsave("results/img/pca13_roots_png.png", pca13_roots,
       width = 20, height = 20, units = "cm", dpi = 300)

### PC2 & PC3 ----
pca23_roots <- pca_sites_roots |> 
  ggplot(aes(x = PC2, y = PC3, 
             colour = species)) +
  geom_point(aes(shape = factor(elevation_m_asl)), size = 4) +
  scale_shape_manual(values = elevation_shapes, name = "Elevation (m asl)") +  # Adjust the number of shapes to match your elevation count
  stat_ellipse(aes(group = species, colour = species), size = 0.8) +
  geom_segment(data = pca_traits_roots,
               aes(x = 0, y = 0, xend = PC2, yend = PC3),
               arrow = arrow(length = unit(0.5, "cm")),
               size = 1,
               colour = "grey20",
               inherit.aes = FALSE) +
  geom_text_repel(data = pca_traits_roots,
                  aes(x = PC2 * 1.1, y = PC3 * 1.1, label = traits),
                  size = 4,
                  fontface = "bold",
                  inherit.aes = FALSE, 
                  colour = "black") +
  coord_equal() +
  scale_colour_manual(values = species_colors, name = "Species",
                      labels = species_labels) +
  labs(x = glue("PCA2 ({round(e_B_roots[2] * 100, 1)}%)"),
       y = glue("PCA3 ({round(e_B_roots[3] * 100, 1)}%)")) +
  theme_bw(base_size = 16) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.ticks = element_line(linewidth = 1),
    axis.line = element_line(linewidth = 1, colour = "black"),
    legend.text = element_text(size = 16),
    plot.margin = margin(2,2,2,2),
    aspect.ratio = 1) +
  coord_cartesian(xlim = xlim_equal_roots)

pca23_roots

#### Save plots ----
# ggsave("results/img/pca23_roots_tiff.tiff", pca23_roots,
#        width = 20, height = 20, units = "cm", dpi = 300)
ggsave("results/img/pca23_roots_png.png", pca23_roots,
       width = 20, height = 20, units = "cm", dpi = 300)


plots_pca_roots <- c(pca12_roots, pca13_roots, pca23_roots)
pca_roots_legend <- get_legend(plots_pca_roots[[1]] + theme(legend.position = "right"))
plots_pca_roots_nolegend <- lapply(plots_pca_roots, function(p) p + theme(legend.position = "none"))


plot_roots_pca <- plot_grid(plots_pca_roots_nolegend[[1]], plots_pca_roots_nolegend[[2]],
                          plots_pca_roots_nolegend[[3]], pca_roots_legend,
                          ncol = 2,
                          rel_widths = c(1, 1), 
                          rel_heights = c(1, 1),
                          labels = c("A.", "B.", "C.", ""),
                          align = "hv",
                          axis = "tblr") + 
  theme(plot.background = element_rect(fill = "white", colour = NA))

plot_roots_pca 

#### Save all plots ----
ggsave("results/img/plot_roots_pca_tiff.tiff", plot_roots_pca,
       width = 25, height = 25, units = "cm", dpi = 300)
ggsave("results/img/plot_roots_pca_png.png", plot_roots_pca,
       width = 25, height = 25, units = "cm", dpi = 300)

# 5. PCA Leaf ----
## PCA using vegan (rda with no constraints = PCA) ----
## Remember to scale to unit variance
pca_leaf <- trait_data_wide |> 
  select(sla, ldmc, leaf_thickness) |> 
  rename(
    LT = leaf_thickness, 
    LDMC = ldmc,
    SLA = sla) |> 
  rda(scale = TRUE)
summary(pca_leaf)

## How many PCs should we keep? ----
### Extract eigenvalues and their percentage
ev_leaf <- eigenvals(pca_leaf)
e_B_leaf <- eigenvals(pca_leaf)/sum(eigenvals(pca_leaf))
### Extract PCs that explain 90% of the variation in the data
k_leaf <- which(cumsum(ev_leaf) / sum(ev_leaf) >= 0.9)[1]

## Tables for ind. loadings and traits loadings ----
## Create a table containing the loading of each invidual (sites) in each PC

pca_sites_leaf <- as_tibble(bind_cols(
  trait_data_wide |> 
    mutate(species = fct_relevel(species, "Themeda triandra", after = 2)) |> 
    select(id, species, family, growth_form, elevation_m_asl),
  scores(pca_leaf, display = "sites", choices = 1:k_leaf, scaling = 2)))

pca_traits_leaf <- scores(pca_leaf, display = "species", choices = 1:k_leaf, scaling = 2) |> 
  as.data.frame() |> 
  rownames_to_column(var = "traits") 

## PCA significance
PCAsignificance(pca_leaf)
plot1_leaf <- ordiplot(pca_leaf, choices=c(1,2), scaling=1)
ordiequilibriumcircle(pca_leaf,plot1_leaf) #Which traits are more important in each PC

## PERMANOVA  ----
### Matrix of distances ----
leaf_distance <- dist(pca_sites_leaf |> 
                         select(starts_with("PC")))
### Species alone ----
set.seed(1)
adonis2(leaf_distance ~ species, data = pca_sites_leaf, permutations = 4999)

# Multivarite homogeneity
leaf_sp_bd <- betadisper(leaf_distance, pca_sites_leaf$species)
anova(leaf_sp_bd)
permutest(leaf_sp_bd, 999)
plot(leaf_sp_bd)

# Plot showing differences among spp
score_leaf_spp <- scores(pca_leaf, display = "sites", scaling = 1, choices = 1:2)
plot(score_leaf_spp, type = "n") 
points(score_leaf_spp, col = as.integer(pca_sites_leaf$species), pch = 19)
ordiellipse(score_leaf_spp, pca_sites_leaf$species, kind = "se", conf = 0.95, draw = "polygon",
            col = 1:5, border = 1:5, label = TRUE)

#### Pairwise comparisons among species ----
set.seed(1)
pairwise.adonis2(leaf_distance ~ species, 
                 data = as.data.frame(pca_sites_leaf),
                 permutations = 4999, 
                 p.adjust.m = "holm")

### Elevation alone ----
set.seed(1)
adonis2(leaf_distance ~ elevation_m_asl, data = pca_sites_leaf, permutations = 4999)

# Multivarite homogeneity
leaf_ele_bd <- betadisper(leaf_distance, pca_sites_leaf$elevation_m_asl)
perm_leaf_ele_bd  <- permutest(leaf_ele_bd, permutations = 4999)
anova(leaf_ele_bd)
permutest(leaf_ele_bd, 999)
plot(leaf_ele_bd)

# Plot showing differences among elevations
score_leaf_ele <- scores(pca_leaf, display = "sites", scaling = 1, choices = 1:2)
plot(score_leaf_ele, type = "n") 
points(score_leaf_ele, col = as.integer(pca_sites_leaf$elevation_m_asl), pch = 19)
ordiellipse(score_leaf_ele, pca_sites_leaf$elevation_m_asl, kind = "se", conf = 0.95, draw = "polygon",
            col = 1:5, border = 1:5, label = TRUE)

#### Pairwise comparisons among elevations ----
set.seed(1)
pairwise.adonis2(leaf_distance ~ elevation_m_asl, 
                 data = as.data.frame(pca_sites_leaf),
                 permutations = 4999, 
                 p.adjust.m = "holm")

### Marginal differences species + elevation ----
set.seed(1)
adonis2(leaf_distance ~ species + elevation_m_asl, data = as.data.frame(pca_sites_leaf), 
        permutations = 4999,
        by = "margin")


pairwise.adonis2(leaf_distance ~ species + elevation_m_asl, 
                 data = as.data.frame(pca_sites_leaf),
                 p.adjust.m = "holm")


### Species pooled “within-elevation” test ----
set.seed(1)
adonis2(leaf_distance ~ species, data = as.data.frame(pca_sites_leaf), permutations = 4999,
        strata = pca_sites_leaf$elevation_m_asl)

#### Pairwise comparisons among species ----
set.seed(1)
pairwise.adonis2(leaf_distance ~ species, 
                 data = as.data.frame(pca_sites_leaf),
                 strata = "elevation_m_asl", 
                 p.adjust.m = "holm")


### Elevation pooled “within-species” test ----
set.seed(1)
adonis2(leaf_distance ~ elevation_m_asl, data = as.data.frame(pca_sites_leaf), 
        permutations = 4999,
        strata = pca_sites_leaf$species)

#### Pairwise comparisons among elevations within species ----
set.seed(1)
pairwise.adonis2(leaf_distance ~ elevation_m_asl, 
                 data = as.data.frame(pca_sites_leaf),
                 strata = "species", 
                 p.adjust.m = "holm")

## PCA Plot ----
xlim_equal_leaf <- c(pca_traits_leaf$PC1, pca_traits_leaf$PC2) |> 
  abs() |> 
  max(na.rm = TRUE) |> 
  (\(m) c(-m, m))()

### PC1 & PC2 ----
pca12_leaf <- pca_sites_leaf |> 
  ggplot(aes(x = PC1, y = PC2, 
             colour = species)) +
  geom_point(aes(shape = factor(elevation_m_asl)), size = 4) +
  scale_shape_manual(values = elevation_shapes, name = "Elevation (m asl)") +  # Adjust the number of shapes to match your elevation count
  stat_ellipse(aes(group = species, colour = species), size = 0.8) +
  geom_segment(data = pca_traits_leaf,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.5, "cm")),
               size = 1,
               colour = "grey20",
               inherit.aes = FALSE) +
  geom_text_repel(data = pca_traits_leaf,
                  aes(x = PC1 * 1.1, y = PC2 * 1.1, label = traits),
                  size = 4,
                  fontface = "bold",
                  inherit.aes = FALSE, 
                  colour = "black") +
  coord_equal() +
  scale_colour_manual(values = species_colors, name = "Species",
                      labels = species_labels) +
  labs(x = glue("PCA1 ({round(e_B_leaf[1] * 100, 1)}%)"),
       y = glue("PCA2 ({round(e_B_leaf[2] * 100, 1)}%)")) +
  theme_bw(base_size = 16) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.ticks = element_line(linewidth = 1),
    axis.line = element_line(linewidth = 1, colour = "black"),
    legend.text = element_text(size = 16),
    plot.margin = margin(2,2,2,2),
    aspect.ratio = 1) +
  coord_cartesian(xlim = xlim_equal_leaf)

pca12_leaf

#### Save plots ----
# ggsave("results/img/pca12_leaf_tiff.tiff", pca12_leaf,
#        width = 20, height = 20, units = "cm", dpi = 300)
ggsave("results/img/pca12_leaf_png.png", pca12_leaf,
       width = 20, height = 20, units = "cm", dpi = 300)

# 6. PCA plant size ----
## PCA using vegan (rda with no constraints = PCA) ----
## Remember to scale to unit variance
pca_plants <- trait_data_wide |>
  select(root_depth, veg_height, bgb_agb) |>
  rename(
    VHeight = veg_height,
    RDepth = root_depth,
    `BG:AG` = bgb_agb) |>
  rda(scale = TRUE)
summary(pca_plants)

## How many PCs should we keep? ----
### Extract eigenvalues and their percentage
ev_plants <- eigenvals(pca_plants)
e_B_plants <- eigenvals(pca_plants)/sum(eigenvals(pca_plants))
### Extract PCs that explain 90% of the variation in the data
k_plants <- which(cumsum(ev_plants) / sum(ev_plants) >= 0.9)[1]

## Tables for ind. loadings and traits loadings ----
## Create a table containing the loading of each invidual (sites) in each PC

pca_sites_plants <- as_tibble(bind_cols(
  trait_data_wide |>
    mutate(species = fct_relevel(species, "Themeda triandra", after = 2)) |>
    select(id, species, family, growth_form, elevation_m_asl),
  scores(pca_plants, display = "sites", choices = 1:k_plants, scaling = 2)))

pca_traits_plants <- scores(pca_plants, display = "species", choices = 1:k_plants, scaling = 2) |>
  as.data.frame() |>
  rownames_to_column(var = "traits")

## PCA significance
PCAsignificance(pca_plants)
plot1_plants <- ordiplot(pca_plants, choices=c(1,2), scaling=1)
ordiequilibriumcircle(pca_plants,plot1_plants) #Which traits are more important in each PC

## PERMANOVA  ----
### Matrix of distances ----
plants_distance <- dist(pca_sites_plants |>
                        select(starts_with("PC")))
### Species alone ----
set.seed(1)
adonis2(plants_distance ~ species, data = pca_sites_plants, permutations = 4999)

# Multivarite homogeneity
plants_sp_bd <- betadisper(plants_distance, pca_sites_plants$species)
anova(plants_sp_bd)
permutest(plants_sp_bd, 999)
plot(plants_sp_bd)

# Plot showing differences among spp
scores_plants_spp <- scores(pca_plants, display = "sites", scaling = 1, choices = 1:2)
plot(scores_plants_spp, type = "n")
points(scores_plants_spp, col = as.integer(pca_sites_plants$species), pch = 19)
ordiellipse(scores_plants_spp, pca_sites_plants$species, kind = "se", conf = 0.95, draw = "polygon",
            col = 1:5, border = 1:5, label = TRUE)

#### Pairwise comparisons among species ----
set.seed(1)
pairwise.adonis2(plants_distance ~ species,
                 data = as.data.frame(pca_sites_plants),
                 permutations = 4999,
                 p.adjust.m = "holm")

### Elevation alone ----
set.seed(1)
adonis2(plants_distance ~ elevation_m_asl, data = pca_sites_plants, permutations = 4999)

# Multivarite homogeneity
plants_ele_bd <- betadisper(plants_distance, pca_sites_plants$elevation_m_asl)
perm_plants_ele_bd  <- permutest(plants_ele_bd, permutations = 4999)
anova(plants_ele_bd)
permutest(plants_ele_bd, 999)
plot(plants_ele_bd)

# Plot showing differences among elevations
scores_plants_ele <- scores(pca_plants, display = "sites", scaling = 1, choices = 1:2)
plot(scores_plants_ele, type = "n")
points(scores_plants_ele, col = as.integer(pca_sites_plants$elevation_m_asl), pch = 19)
ordiellipse(scores_plants_ele, pca_sites_plants$elevation_m_asl, kind = "se", conf = 0.95, draw = "polygon",
            col = 1:5, border = 1:5, label = TRUE)

#### Pairwise comparisons among elevations ----
set.seed(1)
pairwise.adonis2(plants_distance ~ elevation_m_asl,
                 data = as.data.frame(pca_sites_plants),
                 permutations = 4999,
                 p.adjust.m = "holm")

### Marginal differences species + elevation ----
set.seed(1)
adonis2(plants_distance ~ species + elevation_m_asl, data = as.data.frame(pca_sites_plants),
        permutations = 4999,
        by = "margin")


pairwise.adonis2(plants_distance ~ species + elevation_m_asl,
                 data = as.data.frame(pca_sites_plants),
                 p.adjust.m = "holm")


### Species pooled “within-elevation” test ----
set.seed(1)
adonis2(plants_distance ~ species, data = as.data.frame(pca_sites_plants), permutations = 4999,
        strata = pca_sites_plants$elevation_m_asl)

#### Pairwise comparisons among species ----
set.seed(1)
pairwise.adonis2(plants_distance ~ species,
                 data = as.data.frame(pca_sites_plants),
                 strata = "elevation_m_asl",
                 p.adjust.m = "holm")


### Elevation pooled “within-species” test ----
set.seed(1)
adonis2(plants_distance ~ elevation_m_asl, data = as.data.frame(pca_sites_plants),
        permutations = 4999,
        strata = pca_sites_plants$species)

#### Pairwise comparisons among elevations within species ----
set.seed(1)
pairwise.adonis2(plants_distance ~ elevation_m_asl,
                 data = as.data.frame(pca_sites_plants),
                 strata = "species",
                 p.adjust.m = "holm")

## PCA Plot ----
## Color for species
xlim_equal_plants <- c(pca_traits_plants$PC2, pca_traits_plants$PC3) |>
  abs() |>
  max(na.rm = TRUE) |>
  (\(m) c(-m, m))()

### PC1 & PC2 ----
pca12_plants <- pca_sites_plants |>
  ggplot(aes(x = PC1, y = PC2,
             colour = species)) +
  geom_point(aes(shape = factor(elevation_m_asl)), size = 4) +
  scale_shape_manual(values = elevation_shapes, name = "Elevation (m asl)") +  # Adjust the number of shapes to match your elevation count
  stat_ellipse(aes(group = species, colour = species), size = 0.8) +
  geom_segment(data = pca_traits_plants,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.5, "cm")),
               size = 1,
               colour = "grey20",
               inherit.aes = FALSE) +
  geom_text_repel(data = pca_traits_plants,
                  aes(x = PC1 * 1.1, y = PC2 * 1.1, label = traits),
                  size = 4,
                  fontface = "bold",
                  inherit.aes = FALSE,
                  colour = "black") +
  coord_equal() +
  scale_colour_manual(values = species_colors, name = "Species",
                      labels = species_labels) +
  labs(x = glue("PCA1 ({round(e_B_plants[1] * 100, 1)}%)"),
       y = glue("PCA2 ({round(e_B_plants[2] * 100, 1)}%)")) +
  theme_bw(base_size = 16) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.ticks = element_line(linewidth = 1),
    axis.line = element_line(linewidth = 1, colour = "black"),
    legend.text = element_text(size = 16),
    plot.margin = margin(2,2,2,2),
    aspect.ratio = 1) +
  coord_cartesian(xlim = xlim_equal_plants)

pca12_plants

#### Save plots ----
# ggsave("results/img/pca12_plants_tiff.tiff", pca12_plants,
#        width = 20, height = 20, units = "cm", dpi = 300)
ggsave("results/img/pca12_plants_png.png", pca12_plants,
       width = 20, height = 20, units = "cm", dpi = 300)

### PC1 & PC3 ----
pca13_plants <- pca_sites_plants |>
  ggplot(aes(x = PC1, y = PC3,
             colour = species)) +
  geom_point(aes(shape = factor(elevation_m_asl)), size = 4) +
  scale_shape_manual(values = elevation_shapes, name = "Elevation (m asl)") +  # Adjust the number of shapes to match your elevation count
  stat_ellipse(aes(group = species, colour = species), size = 0.8) +
  geom_segment(data = pca_traits_plants,
               aes(x = 0, y = 0, xend = PC1, yend = PC3),
               arrow = arrow(length = unit(0.5, "cm")),
               size = 1,
               colour = "grey20",
               inherit.aes = FALSE) +
  geom_text_repel(data = pca_traits_plants,
                  aes(x = PC1 * 1.1, y = PC3 * 1.1, label = traits),
                  size = 4,
                  fontface = "bold",
                  inherit.aes = FALSE,
                  colour = "black") +
  coord_equal() +
  scale_colour_manual(values = species_colors, name = "Species",
                      labels = species_labels) +
  labs(x = glue("PCA1 ({round(e_B_plants[1] * 100, 1)}%)"),
       y = glue("PCA3 ({round(e_B_plants[3] * 100, 1)}%)")) +
  theme_bw(base_size = 16) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.ticks = element_line(linewidth = 1),
    axis.line = element_line(linewidth = 1, colour = "black"),
    legend.text = element_text(size = 16),
    plot.margin = margin(2,2,2,2),
    aspect.ratio = 1) +
  coord_cartesian(xlim = xlim_equal_plants)

pca13_plants

#### Save plots ----
# ggsave("results/img/pca13_plants_tiff.tiff", pca13_plants,
#        width = 20, height = 20, units = "cm", dpi = 300)
ggsave("results/img/pca13_plants_png.png", pca13_plants,
       width = 20, height = 20, units = "cm", dpi = 300)

### PC2 & PC3 ----
pca23_plants <- pca_sites_plants |>
  ggplot(aes(x = PC2, y = PC3,
             colour = species)) +
  geom_point(aes(shape = factor(elevation_m_asl)), size = 4) +
  scale_shape_manual(values = elevation_shapes, name = "Elevation (m asl)") +  # Adjust the number of shapes to match your elevation count
  stat_ellipse(aes(group = species, colour = species), size = 0.8) +
  geom_segment(data = pca_traits_plants,
               aes(x = 0, y = 0, xend = PC2, yend = PC3),
               arrow = arrow(length = unit(0.5, "cm")),
               size = 1,
               colour = "grey20",
               inherit.aes = FALSE) +
  geom_text_repel(data = pca_traits_plants,
                  aes(x = PC2 * 1.1, y = PC3 * 1.1, label = traits),
                  size = 4,
                  fontface = "bold",
                  inherit.aes = FALSE,
                  colour = "black") +
  coord_equal() +
  scale_colour_manual(values = species_colors, name = "Species",
                      labels = species_labels) +
  labs(x = glue("PCA2 ({round(e_B_plants[2] * 100, 1)}%)"),
       y = glue("PCA3 ({round(e_B_plants[3] * 100, 1)}%)")) +
  theme_bw(base_size = 16) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.ticks = element_line(linewidth = 1),
    axis.line = element_line(linewidth = 1, colour = "black"),
    legend.text = element_text(size = 16),
    plot.margin = margin(2,2,2,2),
    aspect.ratio = 1) +
  coord_cartesian(xlim = xlim_equal_plants)

pca23_plants


#### Save plots ----
# ggsave("results/img/pca23_plants_tiff.tiff", pca23_plants,
#        width = 20, height = 20, units = "cm", dpi = 300)
ggsave("results/img/pca23_plants_png.png", pca23_plants,
       width = 20, height = 20, units = "cm", dpi = 300)


plots_pca_plants <- c(pca12_plants, pca13_plants, pca23_plants)
pca_plants_legend <- get_legend(plots_pca_plants[[1]] + theme(legend.position = "right"))
plots_pca_plants_nolegend <- lapply(plots_pca_plants, function(p) p + theme(legend.position = "none"))


plot_plants_pca <- plot_grid(plots_pca_plants_nolegend[[1]], plots_pca_plants_nolegend[[2]],
                           plots_pca_plants_nolegend[[3]], pca_plants_legend,
                           ncol = 2,
                           rel_widths = c(1, 1),
                           rel_heights = c(1, 1),
                           labels = c("A.", "B.", "C.", ""),
                           align = "hv",
                           axis = "tblr") +
  theme(plot.background = element_rect(fill = "white", colour = NA))

plot_plants_pca

### Save all plots ----
# ggsave("results/img/plot_plants_pca_tiff.tiff", plot_plants_pca,
#        width = 25, height = 25, units = "cm", dpi = 300)
ggsave("results/img/plot_plants_pca_png.png", plot_plants_pca,
       width = 25, height = 25, units = "cm", dpi = 300)


#### Save PCA outputs for all, roots, leaves, size
write_csv(pca_sites, 'data/output/PCA_all_traits.csv')
write_csv(pca_sites_roots, 'data/output/PCA_root_traits.csv')
write_csv(pca_sites_leaf, 'data/output/PCA_leaf_traits.csv')
write_csv(pca_sites_plants, 'data/output/PCA_size_traits.csv')

