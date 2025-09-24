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

#install.packages("devtools", "tidyverse", "ggplot2", "plotly", "vegan", "ggvegan", "ggrepel", "glue", "viridis") #install if needed
# devtools::install_github("gavinsimpson/ggvegan")
pkgs <- c("devtools", "tidyverse", "ggplot2", "plotly", "vegan", "ggvegan", "ggrepel", "glue", "viridis")
lapply(pkgs, library, character.only = TRUE)
remove(pkgs)

# 2. Load data ----
trait_data <- read_csv("data/processed/v_PFCT7_clean_functional_traits_2023.csv")
trait_data_wide <- read_csv("data/processed/v_PFCT7_clean_functional_traits_2023_wide.csv")

# 3. Prepare data ----
# Select numeric columns for PCA
trait_data_sel <- trait_data_wide |> 
  select("rd","bi","srl","rtd","rdmc","sla","ldmc","bgb_agb")


pca_output <- trait_data_wide |> 
  select(rd, bi, srl, rtd, rdmc, sla, ldmc, bgb_agb) |> 
  rda(scale = TRUE)
summary(pca_output)

pca_sites <- bind_cols(
  trait_data_wide %>%
    select(elevation_m_asl, species),
  fortify(pca_output, display = "sites")
)


pca_traits <- ggvegan::fortify(pca_output, display = "species") |> 
  mutate(Trait = label)

# get eigenvalues
e_B <- eigenvals(pca_output)/sum(eigenvals(pca_output))


pca_sites |> 
  ggplot(aes(x = PC1, y = PC2, 
             colour = species)) +
  geom_point(aes(shape = elevation_m_asl), size = 3) +
  geom_segment(data = pca_traits,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")),
               colour = "grey50",
               inherit.aes = FALSE) +
  geom_text_repel(data = pca_traits,
                  aes(x = PC1 * 1.1, y = PC2 * 1.1, label = Trait),
                  size = 4,
                  inherit.aes = FALSE, colour = "black") +
  coord_equal() +
  stat_ellipse(aes(group = species, colour = species), size = 0.8) +
  scale_colour_viridis_d(end = 0.8, begin = 0.1, option = "A", direction = -1, name = "Species") +
  scale_shape_manual(values = c(16, 1, 2, 3), name = "Elevation (masl)") +  # Adjust the number of shapes to match your elevation count
  labs(x = glue("PCA1 ({round(e_B[1] * 100, 1)}%)"),
       y = glue("PCA2 ({round(e_B[2] * 100, 1)}%)")) +
  theme_bw()
