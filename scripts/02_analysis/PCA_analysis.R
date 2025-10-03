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
          "fixest", "lmtest", "corrplot")
lapply(pkgs, library, character.only = TRUE)
remove(pkgs)

# 2. Load data ----
trait_data <- read_csv("data/processed/v_PFCT7_clean_functional_traits_2023.csv")
trait_data_wide <- read_csv("data/processed/v_PFCT7_clean_functional_traits_2023_wide.csv")


# 3. PCA Plot ----
pca_output <- trait_data_wide |> 
  select(root_depth, veg_height, rd, bi, srl, rtd, rdmc, sla, ldmc, bgb_agb) |> 
  rename(
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

pca_sites <- bind_cols(
  trait_data_wide |> 
    mutate(species = fct_relevel(species, "Themeda triandra", after = 2)) |> 
    select(elevation_m_asl, species), 
  fortify(pca_output, display = "sites"))


pca_traits <- ggvegan::fortify(pca_output, display = "species") |> 
  mutate(Trait = label)

# get eigenvalues
e_B <- eigenvals(pca_output)/sum(eigenvals(pca_output))


species_colors <- c(
  "Eragrostis capensis" = "#42049EFF",
  "Harpochloa falx" = "#8204A7FF",
  "Themeda triandra" = "#B6308BFF",
  "Helichrysum pilosellum" = "#F79143FF",
  "Senecio glaberrimus" = "#FCCE25FF")

pca_sites |> 
  ggplot(aes(x = PC1, y = PC2, 
             colour = species)) +
  geom_point(aes(shape = factor(elevation_m_asl)), size = 4) +
  geom_segment(data = pca_traits,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.5, "cm")),
               size = 1,
               colour = "grey20",
               inherit.aes = FALSE) +
  geom_text_repel(data = pca_traits,
                  aes(x = PC1 * 1.1, y = PC2 * 1.1, label = Trait),
                  size = 4,
                  fontface = "bold",
                  inherit.aes = FALSE, 
                  colour = "black") +
  coord_equal() +
  stat_ellipse(aes(group = species, colour = species), size = 0.8) +
  scale_colour_manual(values = species_colors, name = "Species") +
  scale_shape_manual(values = c(19, 15, 17, 18), name = "Elevation (masl)") +  # Adjust the number of shapes to match your elevation count
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
    legend.text = element_text(size = 16))

