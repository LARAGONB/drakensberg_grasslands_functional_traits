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
# devtools::install_github("gavinsimpson/ggvegan")
pkgs <- c("devtools", "tidyverse", "ggplot2", "plotly", "vegan", "ggvegan", "ggrepel", "glue", "viridis",
          "fixest", "lmtest", "corrplot", "FactoMineR", "factoextra", "BiodiversityR",
          "cowplot")
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
## Create a table containing the loading of each invidual (sites) in each PC

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
plot1 <- ordiplot(pca_output, choices=c(1,2), scaling=1)
ordiequilibriumcircle(pca_output,plot1) #Which traits are more important in each PC


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
## Color for species
species_colors <- c(
  "Eragrostis capensis" = "#42049EFF",
  "Harpochloa falx" = "#8204A7FF",
  "Themeda triandra" = "#B6308BFF",
  "Helichrysum pilosellum" = "#F79143FF",
  "Senecio glaberrimus" = "#FCCE25FF")

xlim_equal <- c(pca_sites$PC2, pca_sites$PC3) |> 
  abs() |> 
  max(na.rm = TRUE) |> 
  (\(m) c(-m, m))()



### PC1 & PC2 ----
pca12_all <- pca_sites |> 
  ggplot(aes(x = PC1, y = PC2, 
             colour = species)) +
  geom_point(aes(shape = factor(elevation_m_asl)), size = 4) +
  scale_shape_manual(values = c(19, 15, 17, 3), name = "Elevation (m asl)") +  # Adjust the number of shapes to match your elevation count
  stat_ellipse(aes(group = species, colour = species), size = 0.8) +
  # stat_ellipse(aes(group = species),
  #              type = "euclid", level = 0.95,
  #              linewidth = 1, show.legend = FALSE) +
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
                      labels = c(
                        "Eragrostis capensis" = "ERCA",
                        "Harpochloa falx" = "HAFA",
                        "Themeda triandra" = "THTR",
                        "Helichrysum pilosellum" = "HEPI",
                        "Senecio glaberrimus" = "SEGL")) +
  labs(x = glue("PCA1 ({round(e_B[1] * 100, 1)}%)"),
       y = glue("PCA2 ({round(e_B[2] * 100, 1)}%)")) +
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.ticks = element_line(linewidth = 1),
    # panel.border = element_blank(), 
    axis.line = element_line(linewidth = 1, colour = "black"),
    # legend.title = element_blank(),
    legend.text = element_text(size = 16),
    plot.margin = margin(2,2,2,2),
    aspect.ratio = 1) +
  coord_cartesian(xlim = xlim_equal)
 
pca12_all

#### Save plots ----
ggsave("results/img/pca12_all_tiff.tiff", pca12_all,
       width = 20, height = 20, units = "cm", dpi = 300)
ggsave("results/img/pca12_all_png.png", pca12_all,
       width = 20, height = 20, units = "cm", dpi = 300)

### PC1 & PC3 ----
pca13_all <- pca_sites |> 
  ggplot(aes(x = PC1, y = PC3, 
             colour = species)) +
  geom_point(aes(shape = factor(elevation_m_asl)), size = 4) +
  scale_shape_manual(values = c(19, 15, 17, 3), name = "Elevation (m asl)") +  # Adjust the number of shapes to match your elevation count
  stat_ellipse(aes(group = species, colour = species), size = 0.8) +
  # stat_ellipse(aes(group = species),
  #              type = "euclid", level = 0.95,
  #              linewidth = 1, show.legend = FALSE) +
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
                      labels = c(
                        "Eragrostis capensis" = "ERCA",
                        "Harpochloa falx" = "HAFA",
                        "Themeda triandra" = "THTR",
                        "Helichrysum pilosellum" = "HEPI",
                        "Senecio glaberrimus" = "SEGL")) +
  labs(x = glue("PCA1 ({round(e_B[1] * 100, 1)}%)"),
       y = glue("PCA3 ({round(e_B[3] * 100, 1)}%)")) +
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.ticks = element_line(linewidth = 1),
    # panel.border = element_blank(), 
    axis.line = element_line(linewidth = 1, colour = "black"),
    # legend.title = element_blank(),
    legend.text = element_text(size = 16),
    plot.margin = margin(2,2,2,2),
    aspect.ratio = 1) +
  coord_cartesian(xlim = xlim_equal)

pca13_all

#### Save plots ----
ggsave("results/img/pca13_all_tiff.tiff", pca13_all,
       width = 20, height = 20, units = "cm", dpi = 300)
ggsave("results/img/pca13_all_png.png", pca13_all,
       width = 20, height = 20, units = "cm", dpi = 300)

### PC2 & PC3 ----
pca23_all <- pca_sites |> 
  ggplot(aes(x = PC2, y = PC3, 
             colour = species)) +
  geom_point(aes(shape = factor(elevation_m_asl)), size = 4) +
  scale_shape_manual(values = c(19, 15, 17, 3), name = "Elevation (m asl)") +  # Adjust the number of shapes to match your elevation count
  stat_ellipse(aes(group = species, colour = species), size = 0.8) +
  # stat_ellipse(aes(group = species),
  #              type = "euclid", level = 0.95,
  #              linewidth = 1, show.legend = FALSE) +
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
                      labels = c(
                        "Eragrostis capensis" = "ERCA",
                        "Harpochloa falx" = "HAFA",
                        "Themeda triandra" = "THTR",
                        "Helichrysum pilosellum" = "HEPI",
                        "Senecio glaberrimus" = "SEGL")) +
  labs(x = glue("PCA2 ({round(e_B[2] * 100, 1)}%)"),
       y = glue("PCA3 ({round(e_B[3] * 100, 1)}%)")) +
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.ticks = element_line(linewidth = 1),
    # panel.border = element_blank(), 
    axis.line = element_line(linewidth = 1, colour = "black"),
    # legend.title = element_blank(),
    legend.text = element_text(size = 16),
    plot.margin = margin(2,2,2,2),
    aspect.ratio = 1) +
  coord_cartesian(xlim = xlim_equal)

pca23_all

#### Save plots ----
ggsave("results/img/pca23_all_tiff.tiff", pca23_all,
       width = 20, height = 20, units = "cm", dpi = 300)
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
pca_roots <- trait_data_wide |> 
  select(rd, bi, srl, rtd, rdmc) |> 
  rename(
    BI = bi,
    RD = rd,
    RDMC = rdmc,
    RTD = rtd,
    SRL = srl) |> 
  rda(scale = TRUE)
summary(pca_roots)

## How many PCs should we keep? ----
### Extract eigenvalues and their percentage
ev_roots <- eigenvals(pca_roots)
e_B_roots <- eigenvals(pca_roots)/sum(eigenvals(pca_roots))
### Extract PCs that explain 90% of the variation in the data
k_roots <- which(cumsum(ev_roots) / sum(ev_roots) >= 0.9)[1]

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
PCAsignificance(pca_roots)
plot1_roots <- ordiplot(pca_roots, choices=c(1,2), scaling=1)
ordiequilibriumcircle(pca_roots,plot1_roots) #Which traits are more important in each PC

## PERMANOVA  ----
### Matrix of distances ----
roots_distance <- dist(pca_sites_roots |> 
                       select(starts_with("PC")))
### Species alone ----
set.seed(1)
adonis2(roots_distance ~ species, data = pca_sites_roots, permutations = 4999)

# Multivarite homogeneity
roots_sp_bd <- betadisper(roots_distance, pca_sites_roots$species)
anova(roots_sp_bd)
permutest(roots_sp_bd, 999)
plot(roots_sp_bd)

# Plot showing differences among spp
score_roots_spp <- scores(pca_roots, display = "sites", scaling = 1, choices = 1:2)
plot(score_roots_spp, type = "n") 
points(score_roots_spp, col = as.integer(pca_sites_roots$species), pch = 19)
ordiellipse(score_roots_spp, pca_sites_roots$species, kind = "se", conf = 0.95, draw = "polygon",
            col = 1:5, border = 1:5, label = TRUE)

#### Pairwise comparisons among species ----
set.seed(1)
pairwise.adonis2(roots_distance ~ species, 
                 data = as.data.frame(pca_sites_roots),
                 permutations = 4999, 
                 p.adjust.m = "holm")

### Elevation alone ----
set.seed(1)
adonis2(roots_distance ~ elevation_m_asl, data = pca_sites_roots, permutations = 4999)

# Multivarite homogeneity
roots_ele_bd <- betadisper(roots_distance, pca_sites_roots$elevation_m_asl)
perm_roots_ele_bd  <- permutest(roots_ele_bd, permutations = 4999)
anova(roots_ele_bd)
permutest(roots_ele_bd, 999)
plot(roots_ele_bd)

# Plot showing differences among elevations
score_roots_ele <- scores(pca_roots, display = "sites", scaling = 1, choices = 1:2)
plot(score_roots_ele, type = "n") 
points(score_roots_ele, col = as.integer(pca_sites_roots$elevation_m_asl), pch = 19)
ordiellipse(score_roots_ele, pca_sites_roots$elevation_m_asl, kind = "se", conf = 0.95, draw = "polygon",
            col = 1:5, border = 1:5, label = TRUE)

#### Pairwise comparisons among elevations ----
set.seed(1)
pairwise.adonis2(roots_distance ~ elevation_m_asl, 
                 data = as.data.frame(pca_sites_roots),
                 permutations = 4999, 
                 p.adjust.m = "holm")

### Marginal differences species + elevation ----
set.seed(1)
adonis2(roots_distance ~ species + elevation_m_asl, data = as.data.frame(pca_sites_roots), 
        permutations = 4999,
        by = "margin")


pairwise.adonis2(roots_distance ~ species + elevation_m_asl, 
                 data = as.data.frame(pca_sites_roots),
                 p.adjust.m = "holm")


### Species pooled “within-elevation” test ----
set.seed(1)
adonis2(roots_distance ~ species, data = as.data.frame(pca_sites_roots), permutations = 4999,
        strata = pca_sites_roots$elevation_m_asl)

#### Pairwise comparisons among species ----
set.seed(1)
pairwise.adonis2(roots_distance ~ species, 
                 data = as.data.frame(pca_sites_roots),
                 strata = "elevation_m_asl", 
                 p.adjust.m = "holm")


### Elevation pooled “within-species” test ----
set.seed(1)
adonis2(roots_distance ~ elevation_m_asl, data = as.data.frame(pca_sites_roots), 
        permutations = 4999,
        strata = pca_sites_roots$species)

#### Pairwise comparisons among elevations within species ----
set.seed(1)
pairwise.adonis2(roots_distance ~ elevation_m_asl, 
                 data = as.data.frame(pca_sites_roots),
                 strata = "species", 
                 p.adjust.m = "holm")

## PCA Plot ----
## Color for species
species_colors <- c(
  "Eragrostis capensis" = "#42049EFF",
  "Harpochloa falx" = "#8204A7FF",
  "Themeda triandra" = "#B6308BFF",
  "Helichrysum pilosellum" = "#F79143FF",
  "Senecio glaberrimus" = "#FCCE25FF")

xlim_equal_roots <- c(pca_traits_roots$PC1, pca_traits_roots$PC2) |> 
  abs() |> 
  max(na.rm = TRUE) |> 
  (\(m) c(-m, m))()



### PC1 & PC2 ----
pca12_roots <- pca_sites_roots |> 
  ggplot(aes(x = PC1, y = PC2, 
             colour = species)) +
  geom_point(aes(shape = factor(elevation_m_asl)), size = 4) +
  scale_shape_manual(values = c(19, 15, 17, 3), name = "Elevation (m asl)") +  # Adjust the number of shapes to match your elevation count
  stat_ellipse(aes(group = species, colour = species), size = 0.8) +
  # stat_ellipse(aes(group = species),
  #              type = "euclid", level = 0.95,
  #              linewidth = 1, show.legend = FALSE) +
  geom_segment(data = pca_traits_roots,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.5, "cm")),
               size = 1,
               colour = "grey20",
               inherit.aes = FALSE) +
  geom_text_repel(data = pca_traits_roots,
                  aes(x = PC1 * 1.1, y = PC2 * 1.1, label = traits),
                  size = 4,
                  fontface = "bold",
                  inherit.aes = FALSE, 
                  colour = "black") +
  coord_equal() +
  scale_colour_manual(values = species_colors, name = "Species",
                      labels = c(
                        "Eragrostis capensis" = "ERCA",
                        "Harpochloa falx" = "HAFA",
                        "Themeda triandra" = "THTR",
                        "Helichrysum pilosellum" = "HEPI",
                        "Senecio glaberrimus" = "SEGL")) +
  labs(x = glue("PCA1 ({round(e_B_roots[1] * 100, 1)}%)"),
       y = glue("PCA2 ({round(e_B_roots[2] * 100, 1)}%)")) +
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.ticks = element_line(linewidth = 1),
    # panel.border = element_blank(), 
    axis.line = element_line(linewidth = 1, colour = "black"),
    # legend.title = element_blank(),
    legend.text = element_text(size = 16),
    plot.margin = margin(2,2,2,2),
    aspect.ratio = 1) +
  coord_cartesian(xlim = xlim_equal_roots)

pca12_roots

#### Save plots ----
ggsave("results/img/pca12_roots_tiff.tiff", pca12_roots,
       width = 20, height = 20, units = "cm", dpi = 300)
ggsave("results/img/pca12_roots_png.png", pca12_roots,
       width = 20, height = 20, units = "cm", dpi = 300)

### PC1 & PC3 ----
pca13_roots <- pca_sites_roots |> 
  ggplot(aes(x = PC1, y = PC3, 
             colour = species)) +
  geom_point(aes(shape = factor(elevation_m_asl)), size = 4) +
  scale_shape_manual(values = c(19, 15, 17, 3), name = "Elevation (m asl)") +  # Adjust the number of shapes to match your elevation count
  stat_ellipse(aes(group = species, colour = species), size = 0.8) +
  # stat_ellipse(aes(group = species),
  #              type = "euclid", level = 0.95,
  #              linewidth = 1, show.legend = FALSE) +
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
                      labels = c(
                        "Eragrostis capensis" = "ERCA",
                        "Harpochloa falx" = "HAFA",
                        "Themeda triandra" = "THTR",
                        "Helichrysum pilosellum" = "HEPI",
                        "Senecio glaberrimus" = "SEGL")) +
  labs(x = glue("PCA1 ({round(e_B_roots[1] * 100, 1)}%)"),
       y = glue("PCA3 ({round(e_B_roots[3] * 100, 1)}%)")) +
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.ticks = element_line(linewidth = 1),
    # panel.border = element_blank(), 
    axis.line = element_line(linewidth = 1, colour = "black"),
    # legend.title = element_blank(),
    legend.text = element_text(size = 16),
    plot.margin = margin(2,2,2,2),
    aspect.ratio = 1) +
  coord_cartesian(xlim = xlim_equal_roots)

pca13_roots

#### Save plots ----
ggsave("results/img/pca13_roots_tiff.tiff", pca13_roots,
       width = 20, height = 20, units = "cm", dpi = 300)
ggsave("results/img/pca13_roots_png.png", pca13_roots,
       width = 20, height = 20, units = "cm", dpi = 300)

### PC2 & PC3 ----
pca23_roots <- pca_sites_roots |> 
  ggplot(aes(x = PC2, y = PC3, 
             colour = species)) +
  geom_point(aes(shape = factor(elevation_m_asl)), size = 4) +
  scale_shape_manual(values = c(19, 15, 17, 3), name = "Elevation (m asl)") +  # Adjust the number of shapes to match your elevation count
  stat_ellipse(aes(group = species, colour = species), size = 0.8) +
  # stat_ellipse(aes(group = species),
  #              type = "euclid", level = 0.95,
  #              linewidth = 1, show.legend = FALSE) +
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
                      labels = c(
                        "Eragrostis capensis" = "ERCA",
                        "Harpochloa falx" = "HAFA",
                        "Themeda triandra" = "THTR",
                        "Helichrysum pilosellum" = "HEPI",
                        "Senecio glaberrimus" = "SEGL")) +
  labs(x = glue("PCA2 ({round(e_B_roots[2] * 100, 1)}%)"),
       y = glue("PCA3 ({round(e_B_roots[3] * 100, 1)}%)")) +
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.ticks = element_line(linewidth = 1),
    # panel.border = element_blank(), 
    axis.line = element_line(linewidth = 1, colour = "black"),
    # legend.title = element_blank(),
    legend.text = element_text(size = 16),
    plot.margin = margin(2,2,2,2),
    aspect.ratio = 1) +
  coord_cartesian(xlim = xlim_equal_roots)

pca23_roots

#### Save plots ----
ggsave("results/img/pca23_roots_tiff.tiff", pca23_roots,
       width = 20, height = 20, units = "cm", dpi = 300)
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
## Color for species
species_colors <- c(
  "Eragrostis capensis" = "#42049EFF",
  "Harpochloa falx" = "#8204A7FF",
  "Themeda triandra" = "#B6308BFF",
  "Helichrysum pilosellum" = "#F79143FF",
  "Senecio glaberrimus" = "#FCCE25FF")

xlim_equal_leaf <- c(pca_traits_leaf$PC1, pca_traits_leaf$PC2) |> 
  abs() |> 
  max(na.rm = TRUE) |> 
  (\(m) c(-m, m))()

### PC1 & PC2 ----
pca12_leaf <- pca_sites_leaf |> 
  ggplot(aes(x = PC1, y = PC2, 
             colour = species)) +
  geom_point(aes(shape = factor(elevation_m_asl)), size = 4) +
  scale_shape_manual(values = c(19, 15, 17, 3), name = "Elevation (m asl)") +  # Adjust the number of shapes to match your elevation count
  stat_ellipse(aes(group = species, colour = species), size = 0.8) +
  # stat_ellipse(aes(group = species),
  #              type = "euclid", level = 0.95,
  #              linewidth = 1, show.legend = FALSE) +
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
                      labels = c(
                        "Eragrostis capensis" = "ERCA",
                        "Harpochloa falx" = "HAFA",
                        "Themeda triandra" = "THTR",
                        "Helichrysum pilosellum" = "HEPI",
                        "Senecio glaberrimus" = "SEGL")) +
  labs(x = glue("PCA1 ({round(e_B_leaf[1] * 100, 1)}%)"),
       y = glue("PCA2 ({round(e_B_leaf[2] * 100, 1)}%)")) +
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.ticks = element_line(linewidth = 1),
    # panel.border = element_blank(), 
    axis.line = element_line(linewidth = 1, colour = "black"),
    # legend.title = element_blank(),
    legend.text = element_text(size = 16),
    plot.margin = margin(2,2,2,2),
    aspect.ratio = 1) +
  coord_cartesian(xlim = xlim_equal_leaf)

pca12_leaf

#### Save plots ----
ggsave("results/img/pca12_leaf_tiff.tiff", pca12_leaf,
       width = 20, height = 20, units = "cm", dpi = 300)
ggsave("results/img/pca12_leaf_png.png", pca12_leaf,
       width = 20, height = 20, units = "cm", dpi = 300)

### PC1 & PC3 ----
# pca13_leaf <- pca_sites_leaf |> 
#   ggplot(aes(x = PC1, y = PC3, 
#              colour = species)) +
#   geom_point(aes(shape = factor(elevation_m_asl)), size = 4) +
#   scale_shape_manual(values = c(19, 15, 17, 3), name = "Elevation (m asl)") +  # Adjust the number of shapes to match your elevation count
#   stat_ellipse(aes(group = species, colour = species), size = 0.8) +
#   #stat_ellipse(aes(group = species),
#   #              type = "euclid", level = 0.95,
#   #              linewidth = 1, show.legend = FALSE) +
#   geom_segment(data = pca_traits_leaf,
#                aes(x = 0, y = 0, xend = PC1, yend = PC3),
#                arrow = arrow(length = unit(0.5, "cm")),
#                size = 1,
#                colour = "grey20",
#                inherit.aes = FALSE) +
#   geom_text_repel(data = pca_traits_leaf,
#                   aes(x = PC1 * 1.1, y = PC3 * 1.1, label = traits),
#                   size = 4,
#                   fontface = "bold",
#                   inherit.aes = FALSE, 
#                   colour = "black") +
#   coord_equal() +
#   scale_colour_manual(values = species_colors, name = "Species",
#                       labels = c(
#                         "Eragrostis capensis" = "ERCA",
#                         "Harpochloa falx" = "HAFA",
#                         "Themeda triandra" = "THTR",
#                         "Helichrysum pilosellum" = "HEPI",
#                         "Senecio glaberrimus" = "SEGL")) +
#   labs(x = glue("PCA1 ({round(e_B_leaf[1] * 100, 1)}%)"),
#        y = glue("PCA3 ({round(e_B_leaf[3] * 100, 1)}%)")) +
#   theme_bw(base_size = 14) +
#   theme(
#     axis.title = element_text(size = 16),
#     axis.text = element_text(size = 16, color = "black"),
#     axis.ticks = element_line(linewidth = 1),
#     # panel.border = element_blank(), 
#     axis.line = element_line(linewidth = 1, colour = "black"),
#     # legend.title = element_blank(),
#     legend.text = element_text(size = 16),
#     plot.margin = margin(2,2,2,2),
#     aspect.ratio = 1) +
#   coord_cartesian(xlim = xlim_equal_leaf)
# 
# pca13_leaf
#
#### Save plots ----
# ggsave("results/img/pca13_leaf_tiff.tiff", pca13_leaf,
#        width = 20, height = 20, units = "cm", dpi = 300)
# ggsave("results/img/pca13_leaf_png.png", pca13_leaf,
#        width = 20, height = 20, units = "cm", dpi = 300)
# 
### PC2 & PC3 ----
# pca23_leaf <- pca_sites_leaf |> 
#   ggplot(aes(x = PC2, y = PC3, 
#              colour = species)) +
#   geom_point(aes(shape = factor(elevation_m_asl)), size = 4) +
#   scale_shape_manual(values = c(19, 15, 17, 3), name = "Elevation (m asl)") +  # Adjust the number of shapes to match your elevation count
#   stat_ellipse(aes(group = species, colour = species), size = 0.8) +
#   # stat_ellipse(aes(group = species),
#   #              type = "euclid", level = 0.95,
#   #              linewidth = 1, show.legend = FALSE) +
#   geom_segment(data = pca_traits_leaf,
#                aes(x = 0, y = 0, xend = PC2, yend = PC3),
#                arrow = arrow(length = unit(0.5, "cm")),
#                size = 1,
#                colour = "grey20",
#                inherit.aes = FALSE) +
#   geom_text_repel(data = pca_traits_leaf,
#                   aes(x = PC2 * 1.1, y = PC3 * 1.1, label = traits),
#                   size = 4,
#                   fontface = "bold",
#                   inherit.aes = FALSE, 
#                   colour = "black") +
#   coord_equal() +
#   scale_colour_manual(values = species_colors, name = "Species",
#                       labels = c(
#                         "Eragrostis capensis" = "ERCA",
#                         "Harpochloa falx" = "HAFA",
#                         "Themeda triandra" = "THTR",
#                         "Helichrysum pilosellum" = "HEPI",
#                         "Senecio glaberrimus" = "SEGL")) +
#   labs(x = glue("PCA2 ({round(e_B_leaf[2] * 100, 1)}%)"),
#        y = glue("PCA3 ({round(e_B_leaf[3] * 100, 1)}%)")) +
#   theme_bw(base_size = 14) +
#   theme(
#     axis.title = element_text(size = 16),
#     axis.text = element_text(size = 16, color = "black"),
#     axis.ticks = element_line(linewidth = 1),
#     # panel.border = element_blank(), 
#     axis.line = element_line(linewidth = 1, colour = "black"),
#     # legend.title = element_blank(),
#     legend.text = element_text(size = 16),
#     plot.margin = margin(2,2,2,2),
#     aspect.ratio = 1) +
#   coord_cartesian(xlim = xlim_equal_leaf)
# 
# pca23_leaf
#
#### Save plots ----
# ggsave("results/img/pca23_leaf_tiff.tiff", pca23_leaf,
#        width = 20, height = 20, units = "cm", dpi = 300)
# ggsave("results/img/pca23_leaf_png.png", pca23_leaf,
#        width = 20, height = 20, units = "cm", dpi = 300)
# 
# 
# plots_pca_leaf <- c(pca12_leaf, pca13_leaf, pca23_leaf)
# pca_leaf_legend <- get_legend(plots_pca_leaf[[1]] + theme(legend.position = "right"))
# plots_pca_leaf_nolegend <- lapply(plots_pca_leaf, function(p) p + theme(legend.position = "none"))
# 
# 
# plot_leaf_pca <- plot_grid(plots_pca_leaf_nolegend[[1]], plots_pca_leaf_nolegend[[2]],
#                             plots_pca_leaf_nolegend[[3]], pca_leaf_legend,
#                             ncol = 2,
#                             rel_widths = c(1, 1), 
#                             rel_heights = c(1, 1),
#                             labels = c("A.", "B.", "C.", ""),
#                             align = "hv",
#                             axis = "tblr") + 
#   theme(plot.background = element_rect(fill = "white", colour = NA))
# 
# plot_leaf_pca 
#
#### Save all plots ----
# ggsave("results/img/plot_leaf_pca_tiff.tiff", plot_leaf_pca,
#        width = 25, height = 25, units = "cm", dpi = 300)
# 
# 
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
species_colors <- c(
  "Eragrostis capensis" = "#42049EFF",
  "Harpochloa falx" = "#8204A7FF",
  "Themeda triandra" = "#B6308BFF",
  "Helichrysum pilosellum" = "#F79143FF",
  "Senecio glaberrimus" = "#FCCE25FF")

xlim_equal_plants <- c(pca_traits_plants$PC2, pca_traits_plants$PC3) |>
  abs() |>
  max(na.rm = TRUE) |>
  (\(m) c(-m, m))()



### PC1 & PC2 ----
pca12_plants <- pca_sites_plants |>
  ggplot(aes(x = PC1, y = PC2,
             colour = species)) +
  geom_point(aes(shape = factor(elevation_m_asl)), size = 4) +
  scale_shape_manual(values = c(19, 15, 17, 3), name = "Elevation (m asl)") +  # Adjust the number of shapes to match your elevation count
  stat_ellipse(aes(group = species, colour = species), size = 0.8) +
  # stat_ellipse(aes(group = species),
  #              type = "euclid", level = 0.95,
  #              linewidth = 1, show.legend = FALSE) +
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
                      labels = c(
                        "Eragrostis capensis" = "ERCA",
                        "Harpochloa falx" = "HAFA",
                        "Themeda triandra" = "THTR",
                        "Helichrysum pilosellum" = "HEPI",
                        "Senecio glaberrimus" = "SEGL")) +
  labs(x = glue("PCA1 ({round(e_B_plants[1] * 100, 1)}%)"),
       y = glue("PCA2 ({round(e_B_plants[2] * 100, 1)}%)")) +
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.ticks = element_line(linewidth = 1),
    # panel.border = element_blank(),
    axis.line = element_line(linewidth = 1, colour = "black"),
    # legend.title = element_blank(),
    legend.text = element_text(size = 16),
    plot.margin = margin(2,2,2,2),
    aspect.ratio = 1) +
  coord_cartesian(xlim = xlim_equal_plants)

pca12_plants

#### Save plots ----
ggsave("results/img/pca12_plants_tiff.tiff", pca12_plants,
       width = 20, height = 20, units = "cm", dpi = 300)
ggsave("results/img/pca12_plants_png.png", pca12_plants,
       width = 20, height = 20, units = "cm", dpi = 300)

### PC1 & PC3 ----
pca13_plants <- pca_sites_plants |>
  ggplot(aes(x = PC1, y = PC3,
             colour = species)) +
  geom_point(aes(shape = factor(elevation_m_asl)), size = 4) +
  scale_shape_manual(values = c(19, 15, 17, 3), name = "Elevation (m asl)") +  # Adjust the number of shapes to match your elevation count
  stat_ellipse(aes(group = species, colour = species), size = 0.8) +
  # stat_ellipse(aes(group = species),
  #              type = "euclid", level = 0.95,
  #              linewidth = 1, show.legend = FALSE) +
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
                      labels = c(
                        "Eragrostis capensis" = "ERCA",
                        "Harpochloa falx" = "HAFA",
                        "Themeda triandra" = "THTR",
                        "Helichrysum pilosellum" = "HEPI",
                        "Senecio glaberrimus" = "SEGL")) +
  labs(x = glue("PCA1 ({round(e_B_plants[1] * 100, 1)}%)"),
       y = glue("PCA3 ({round(e_B_plants[3] * 100, 1)}%)")) +
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.ticks = element_line(linewidth = 1),
    # panel.border = element_blank(),
    axis.line = element_line(linewidth = 1, colour = "black"),
    # legend.title = element_blank(),
    legend.text = element_text(size = 16),
    plot.margin = margin(2,2,2,2),
    aspect.ratio = 1) +
  coord_cartesian(xlim = xlim_equal_plants)

pca13_plants

#### Save plots ----
ggsave("results/img/pca13_plants_tiff.tiff", pca13_plants,
       width = 20, height = 20, units = "cm", dpi = 300)
ggsave("results/img/pca13_plants_png.png", pca13_plants,
       width = 20, height = 20, units = "cm", dpi = 300)

### PC2 & PC3 ----
pca23_plants <- pca_sites_plants |>
  ggplot(aes(x = PC2, y = PC3,
             colour = species)) +
  geom_point(aes(shape = factor(elevation_m_asl)), size = 4) +
  scale_shape_manual(values = c(19, 15, 17, 3), name = "Elevation (m asl)") +  # Adjust the number of shapes to match your elevation count
  stat_ellipse(aes(group = species, colour = species), size = 0.8) +
  # stat_ellipse(aes(group = species),
  #              type = "euclid", level = 0.95,
  #              linewidth = 1, show.legend = FALSE) +
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
                      labels = c(
                        "Eragrostis capensis" = "ERCA",
                        "Harpochloa falx" = "HAFA",
                        "Themeda triandra" = "THTR",
                        "Helichrysum pilosellum" = "HEPI",
                        "Senecio glaberrimus" = "SEGL")) +
  labs(x = glue("PCA2 ({round(e_B_plants[2] * 100, 1)}%)"),
       y = glue("PCA3 ({round(e_B_plants[3] * 100, 1)}%)")) +
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.ticks = element_line(linewidth = 1),
    # panel.border = element_blank(),
    axis.line = element_line(linewidth = 1, colour = "black"),
    # legend.title = element_blank(),
    legend.text = element_text(size = 16),
    plot.margin = margin(2,2,2,2),
    aspect.ratio = 1) +
  coord_cartesian(xlim = xlim_equal_plants)

pca23_plants


#### Save plots ----
ggsave("results/img/pca23_plants_tiff.tiff", pca23_plants,
       width = 20, height = 20, units = "cm", dpi = 300)
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
ggsave("results/img/plot_plants_pca_tiff.tiff", plot_plants_pca,
       width = 25, height = 25, units = "cm", dpi = 300)
ggsave("results/img/plot_plants_pca_png.png", plot_plants_pca,
       width = 25, height = 25, units = "cm", dpi = 300)


# 7. MFA Analysis ----

data_ordered <- trait_data_wide |> 
  select(leaf_thickness, ldmc, sla, 
         bi, rd, rdmc, rtd, srl, 
         veg_height, root_depth, bgb_agb) |> 
  rename(LT = leaf_thickness, 
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

grp_sizes <- c(3, 5, 3) #3 leaf traits, 5 root traits, 3 plant size traits
ncp_max <- min(nrow(data_ordered) - 1L, sum(c(3,5,3)))

## MFA ----
mfa_all <- FactoMineR::MFA(
  data_ordered,
  group = grp_sizes,
  type  = c("s", "s", "s"),                 # per-group types
  name.group = c("Leaf", "Roots", "Plant size"),
  ncp           = ncp_max,
  graph         = FALSE
)

### Screeplot ----
fviz_screeplot(mfa_all)

### Fviz contributions Dim 1 ----
fviz_contrib(mfa_all, choice = "quanti.var", axes = 1, top = 20,
             palette = "jco")

### Fviz contributions Dim 2 ----
fviz_contrib(mfa_all, choice = "quanti.var", axes = 2, top = 20,
             palette = "jco")

### Fviz MFA plot ----
fviz_mfa_var(mfa_all, "quanti.var", palette = "jco", 
             col.var.sup = "violet", repel = TRUE)


### MFA Contribution table ----
all_contrib <- enframe(mfa_all$quanti.var$contrib,
                       name = "variable") |> 
  deframe() |> 
  as.data.frame() |> 
  rownames_to_column(var = "traits") |> 
  left_join(x = _, 
            y = mfa_all$summary.quanti |> 
              select(group, variable),
            by = join_by(traits == variable)) |> 
  mutate(
    group = case_when(
      group == 1 ~ "Leaf",
      group == 2 ~ "Roots",
      group == 3 ~ "Plant size")) |> 
  relocate(group, .after = "traits")

### Contri Dim 1 plot ----
all_contrib_dim1 <- all_contrib |> 
  mutate(traits = fct_reorder(traits, Dim.1, .desc = TRUE)) |> 
  mutate(group = fct_relevel(group, "Leaf", "Roots", "Plant size"))|> 
  ggplot(aes(x = traits, y = Dim.1)) +
  geom_col(aes(fill = group)) +
  scale_fill_manual(values = c(Leaf = "#117733", Roots = "#7F3B08", "Plant size" = "#0072B2")) +
  geom_hline(yintercept = mean(all_contrib$Dim.1), linetype = "dashed") +
  labs(y = "Contributions to Dim-1 (%)",
       x = "",
       fill = "") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.position = c(1, 1),
    legend.justification.inside = c(1,1))

all_contrib_dim1

#### Save plots ----
ggsave("results/img/all_contrib_dim1_tiff.tiff", all_contrib_dim1,
       width = 20, height = 10, units = "cm", dpi = 300)
ggsave("results/img/all_contrib_dim1_png.png", all_contrib_dim1,
       width = 20, height = 10, units = "cm", dpi = 300)

### Contri Dim 2 plots ----
all_contrib_dim2 <- all_contrib |> 
  mutate(traits = fct_reorder(traits, Dim.2, .desc = TRUE)) |> 
  mutate(group = fct_relevel(group, "Leaf", "Roots", "Plant size"))|> 
  ggplot(aes(x = traits, y = Dim.2)) +
  geom_col(aes(fill = group)) +
  scale_fill_manual(values = c(Leaf = "#117733", Roots = "#7F3B08", "Plant size" = "#0072B2")) +
  geom_hline(yintercept = mean(all_contrib$Dim.1), linetype = "dashed") +
  labs(y = "Contributions to Dim-2 (%)",
       x = "",
       fill = "") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.position = c(1, 1),
    legend.justification.inside = c(1,1))

all_contrib_dim2


#### Save plots ----
ggsave("results/img/all_contrib_dim2_tiff.tiff", all_contrib_dim2,
       width = 20, height = 10, units = "cm", dpi = 300)
ggsave("results/img/all_contrib_dim2_png.png", all_contrib_dim2,
       width = 20, height = 10, units = "cm", dpi = 300)

#### Save both contrib plots ----

all_contrib_12 <- plot_grid(all_contrib_dim1, all_contrib_dim2  + theme(legend.position = "none"),
                            nrow = 2,
                            labels = c('A.',"B."),
                            align = "hv")

all_contrib_12

ggsave("results/img/all_contrib_12_tiff.tiff", all_contrib_12,
       width = 20, height = 22, units = "cm", dpi = 300)
ggsave("results/img/all_contrib_12_png.png", all_contrib_12,
       width = 20, height = 22, units = "cm", dpi = 300)


# A. Extract ind. coordinate from MFA ----
ind_coord <- as_tibble(mfa_all$ind$coord)
ind_coord$Row <- row.names(ind_coord)

# B. Create meta data table
meta_data <- trait_data_wide |>
  mutate(Row = as.character(dplyr::row_number())) |> 
  select(Row, id, species, family, growth_form, elevation_m_asl)

# C. Join tables 

mfa_table <- meta_data |> 
  left_join(ind_coord, by = "Row")

# D. Variable loadings (traits arrows)
var <- as_tibble(rownames_to_column(as.data.frame(mfa_all$quanti.var$coord),
                                    var = "traits"))

# E. Axis labels
ev <- mfa_all$eig[,2]
ax1 <- sprintf("Dim 1 (%.1f%%)", ev[1])
ax2 <- sprintf("Dim 2 (%.1f%%)", ev[2])

# F. Scale arrows to fit panel
rng_ind1 <- range(mfa_table$Dim.1)
rng_ind2 <- range(mfa_table$Dim.2)

rng_var1 <- range(var$Dim.1)
rng_var2 <- range(var$Dim.2)

sf <- 0.9 * min(diff(rng_ind1)/diff(rng_var1),
                diff(rng_ind2)/diff(rng_var2))

var_scaled <- transform(var, 
                        xend = Dim.1 * sf, 
                        yend = Dim.2 * sf)
# D. Plot ----

ggplot(mfa_table, aes(x = Dim.1, y = Dim.2, 
                      color = species)) +
  geom_point(aes(shape = as.factor(elevation_m_asl)), size = 3, alpha = 0.9) +
  stat_ellipse(aes(group = species),
               type = "norm", level = 0.95, 
               linewidth = 1, show.legend = FALSE) +
  geom_segment(data = var_scaled,
               aes(x = 0, y = 0, xend = xend, yend = yend),
               inherit.aes = FALSE,
               arrow = arrow(length = unit(0.18, "cm")),
               color = "gray30", linewidth = 0.6) +
  ggrepel::geom_text_repel(
    data = var_scaled, inherit.aes = FALSE,
    aes(x = xend, y = yend, label = traits),
    color = "gray20", size = 3) +
  coord_equal() +
  scale_colour_manual(values = species_colors, name = "Species") +
  scale_shape_manual(values = c(19, 15, 17, 3), name = "Elevation (masl)") +
  labs(x = ax1,
       y = ax2)
