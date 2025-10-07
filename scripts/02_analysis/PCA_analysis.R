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
          "fixest", "lmtest", "corrplot", "FactoMineR", "factoextra")
lapply(pkgs, library, character.only = TRUE)
remove(pkgs)

# 2. Load data ----
trait_data <- read_csv("data/processed/v_PFCT7_clean_functional_traits_2023.csv")
trait_data_wide <- read_csv("data/processed/v_PFCT7_clean_functional_traits_2023_wide.csv")

# 3. PCA Plot ----
pca_output <- trait_data_wide |> 
  select(root_depth, veg_height, rd, bi, srl, rtd, rdmc, sla, ldmc, leaf_thickness, bgb_agb) |> 
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
  # stat_ellipse(aes(group = species),
  #              type = "norm", level = 0.95,
  #              linewidth = 1, show.legend = FALSE) +
  scale_colour_manual(values = species_colors, name = "Species",
                      labels = c(
                        "Eragrostis capensis" = "ERCA",
                        "Harpochloa falx" = "HAFA",
                        "Themeda triandra" = "THTR",
                        "Helichrysum pilosellum" = "HEPI",
                        "Senecio glaberrimus" = "SEGL")) +
  scale_shape_manual(values = c(19, 15, 17, 3), name = "Elevation (m asl)") +  # Adjust the number of shapes to match your elevation count
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
 
PCAsignificance(pca_output)
plot1 <- ordiplot(pca_output, choices=c(1,2), scaling=1)
ordiequilibriumcircle(pca_output,plot1)

vegan::
anova.cca(pca_output, step = 1000)




# 4. MFA Analysis ----

data_ordered <- trait_data_wide |> 
  select(leaf_thickness, ldmc, sla, 
         bi, rd, rdmc, rtd, srl, 
         veg_height, root_depth, bgb_agb) |> 
  rename(LT = leaf_thickness, LDMC = ldmc, SLA = sla, BI = bi, RD = rd, RDMC = rdmc, RTD = rtd, 
         SRL = srl, VHeight = veg_height, RDepth = root_depth, `BG:AG` = bgb_agb)

grp_sizes <- c(3, 5, 3) #3 leaf traits, 5 root traits, 3 plant size traits
ncp_max <- min(nrow(data_ordered) - 1L, sum(c(3,5,3)))
 
mfa_all <- FactoMineR::MFA(
  data_ordered,
  group = grp_sizes,
  type  = c("s", "s", "s"),                 # per-group types
  name.group = c("Leaf", "Roots", "Plant size"),
  ncp           = ncp_max,
  graph         = FALSE
)

#Screeplot
fviz_screeplot(mfa_all)
# Contribution to the first dimension
fviz_contrib(mfa_all, "group", axes = 1)
# Contribution to the second dimension
fviz_contrib(mfa_all, "group", axes = 2)
# Coordinates
quanti.var <- get_mfa_var(mfa_all, "quanti.var")
head(quanti.var$coord)
# Cos2: quality on the factore map
head(quanti.var$cos2)
# Contributions to the dimensions
head(quanti.var$contrib)
# Contributions to dimension 1
fviz_contrib(mfa_all, choice = "quanti.var", axes = 1, top = 20,
             palette = "jco")
# Contributions to dimension 2
fviz_contrib(mfa_all, choice = "quanti.var", axes = 2, top = 20,
             palette = "jco")

#Plot from fviz
fviz_mfa_var(mfa_all, "quanti.var", palette = "jco", 
             col.var.sup = "violet", repel = TRUE)


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
