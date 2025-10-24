################################################################################
# MFA analyses
################################################################################
#
# Lina Aragón
# linamaragonb@gmail.com
# October 14, 2025
#
# Description
################################################################################

# 1. Load libraries ----

#install.packages("devtools", "tidyverse", "ggplot2", "plotly", "vegan", "ggvegan", "ggrepel", "glue", "viridis", "fixest", "lmtest", "corrplot") #install if needed
# devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
pkgs <- c("devtools", "tidyverse", "ggplot2", "plotly", "vegan", "ggvegan", "ggrepel", "glue", "viridis",
          "fixest", "lmtest", "corrplot", "FactoMineR", "factoextra", "BiodiversityR",
          "cowplot", "pairwiseAdonis", "patchwork", "emmeans", "ggExtra")
lapply(pkgs, library, character.only = TRUE)
remove(pkgs)

# 2. Load data ----
trait_data_wide <- read_csv("data/processed/v_PFCT7_clean_functional_traits_2023_wide.csv")

# 3. MFA for all species ----
## Function ----
process_all <- function(data) {
  
  # Get data for current species (without species column)
  ordered_data <- data |> 
    select(id, species, family, growth_form, elevation_m_asl,  # Keep metadata
           leaf_thickness, ldmc, sla, 
           bi, rd, rdmc, rtd, srl, 
           veg_height, root_depth, bgb_agb) |> 
    rename(LT = leaf_thickness, LDMC = ldmc, SLA = sla, 
           BI = bi, RD = rd, RDMC = rdmc, RTD = rtd, SRL = srl, 
           VHeight = veg_height, RDepth = root_depth, `BG:AG` = bgb_agb)
  
  # Remove missing data but keep track of which rows remain
  all_clean <- na.omit(ordered_data)
  
  if (nrow(all_clean) < 3) {
    cat("  Skipping - insufficient data\n")
    return(NULL)
  }
  
  # Extract metadata for the remaining rows (after na.omit)
  metadata <- all_clean |> 
    select(id, species, family, growth_form, elevation_m_asl)
  
  # Extract only trait data for MFA
  trait_data <- all_clean |> 
    select(LT:last_col())  # From LT to BG:AG
  
  # Calculate ncp
  grp_sizes <- c(3, 5, 3)
  ncp_max_sp <- min(nrow(trait_data) - 1L, sum(grp_sizes) - length(grp_sizes))
  
  # Run MFA on trait data only
  result <- MFA(
    base = trait_data,
    group = grp_sizes,
    type = c("s", "s", "s"),
    name.group = c("Leaf", "Roots", "Plant size"),
    ncp = ncp_max_sp,
    graph = FALSE
  )
  
  # Store IDs properly
  rownames(result$ind$coord) <- metadata$id  # Use ID column as rownames
  
  # Create comprehensive coordinates table with all metadata
  result$coord_with_metadata <- result$ind$coord |> 
    as_tibble(rownames = "id") |> 
    left_join(metadata, by = "id") |> 
    relocate(id, species, family, growth_form, elevation_m_asl, .before = everything())
  
  cat("  Success! N =", nrow(trait_data), "\n")
  return(result)
}

### Run the function ----
mfa_all <- process_all(trait_data_wide)

# Access the loadings with metadata
mfa_table <- mfa_all$coord_with_metadata

## MFA Contribution table ----
all_contrib <- mfa_all$quanti.var$contrib |> 
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
  arrange(desc(Dim.1)) |> 
  mutate(traits = factor(traits, levels = traits)) |>
  mutate(group = fct_relevel(group, "Leaf", "Roots", "Plant size"))|> 
  ggplot(aes(x = traits, y = Dim.1)) +
  geom_col(aes(fill = group)) +
  scale_fill_manual(values = group_colors) +
  geom_hline(yintercept = mean(all_contrib$Dim.1), linetype = "dashed") +
  labs(y = "Contributions to Dim-1 (%)",
       x = "",
       fill = "") +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.position = c(1, 1),
    legend.justification.inside = c(1,1))

all_contrib_dim1

### Contri Dim 2 plots ----
all_contrib_dim2 <- all_contrib |> 
  arrange(desc(Dim.2)) |> 
  mutate(traits = factor(traits, levels = traits)) |>
  mutate(group = fct_relevel(group, "Leaf", "Roots", "Plant size"))|> 
  ggplot(aes(x = traits, y = Dim.2)) +
  geom_col(aes(fill = group)) +
  scale_fill_manual(values = group_colors) +
  geom_hline(yintercept = mean(all_contrib$Dim.1), linetype = "dashed") +
  labs(y = "Contributions to Dim-2 (%)",
       x = "",
       fill = "") +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.position = c(1, 1),
    legend.justification.inside = c(1,1))

all_contrib_dim2

#### Save both contrib plots ----

all_contrib_12 <- plot_grid(all_contrib_dim1, all_contrib_dim2  + theme(legend.position = "none"),
                            nrow = 2,
                            labels = c('A.',"B."),
                            align = "hv")

all_contrib_12

# ggsave("results/img/all_contrib_12_tiff.tiff", all_contrib_12,
#        width = 20, height = 22, units = "cm", dpi = 300)
ggsave("results/img/all_contrib_12_png.png", all_contrib_12,
       width = 15, height = 22, units = "cm", dpi = 300)

## MFA Results ----
### Ind. coords table ----
mfa_all_ind <- mfa_table |> 
  mutate(species = fct_relevel(species, "Themeda triandra", after = 2))

write_csv(mfa_all_ind, "data/output/mfa_all_ind.csv")

### Traits loadings (arrows) + scale ----
mfa_all_traits <- {
  
  # Get coordinates
  ind <- mfa_all$ind$coord |> 
    as_tibble()
  var <- mfa_all$quanti.var$coord |> 
    as_tibble(rownames = "traits")
  
  # Get dimension names
  dim_names <- names(ind)
  n_traits <- nrow(var)
  
  # Create scaled dims for common pairs
  scaled_dims <- tibble(.rows = n_traits)
  
  # Dim 1 vs Dim 2 (most common)
  if (length(dim_names) >= 2) {
    rng_ind1 <- range(ind[[dim_names[1]]])
    rng_ind2 <- range(ind[[dim_names[2]]])
    rng_var1 <- range(var[[dim_names[1]]])
    rng_var2 <- range(var[[dim_names[2]]])
    
    sf_12 <- 0.9 * min(diff(rng_ind1)/diff(rng_var1),
                       diff(rng_ind2)/diff(rng_var2))
    
    scaled_dims <- scaled_dims |> 
      bind_cols(
        tibble(
          Dim1_scaled_12 = var[[dim_names[1]]] * sf_12,
          Dim2_scaled_12 = var[[dim_names[2]]] * sf_12
        )
      )
  }
  
  # Dim 1 vs Dim 3 (if exists)
  if (length(dim_names) >= 3) {
    rng_ind1 <- range(ind[[dim_names[1]]])
    rng_ind3 <- range(ind[[dim_names[3]]])
    rng_var1 <- range(var[[dim_names[1]]])
    rng_var3 <- range(var[[dim_names[3]]])
    
    sf_13 <- 0.9 * min(diff(rng_ind1)/diff(rng_var1),
                       diff(rng_ind3)/diff(rng_var3))
    
    scaled_dims <- scaled_dims |> 
      bind_cols(
        tibble(
          Dim1_scaled_13 = var[[dim_names[1]]] * sf_13,
          Dim3_scaled_13 = var[[dim_names[3]]] * sf_13
        )
      )
  }
  
  # Dim 2 vs Dim 3 (if exists)
  if (length(dim_names) >= 3) {
    rng_ind2 <- range(ind[[dim_names[2]]])
    rng_ind3 <- range(ind[[dim_names[3]]])
    rng_var2 <- range(var[[dim_names[2]]])
    rng_var3 <- range(var[[dim_names[3]]])
    
    sf_23 <- 0.9 * min(diff(rng_ind2)/diff(rng_var2),
                       diff(rng_ind3)/diff(rng_var3))
    
    scaled_dims <- scaled_dims |> 
      bind_cols(
        tibble(
          Dim2_scaled_23 = var[[dim_names[2]]] * sf_23,
          Dim3_scaled_23 = var[[dim_names[3]]] * sf_23
        )
      )
  }
  
  # Combine everything 
  scaled_dims |>  
    bind_cols(var) |>  # Add the original variable coordinates
    left_join(
      mfa_all$summary.quanti |> select(group, variable),
      by = join_by(traits == variable)
    ) |> 
    mutate(
      group = case_when(
        group == "1" ~ "Leaf",
        group == "2" ~ "Roots",
        group == "3" ~ "Plant size"),
      group = fct_relevel(group, "Roots", after = 1))|> 
    relocate(traits, group, .before = everything()) 
}

write_csv(mfa_all_traits, "data/output/mfa_all_traits.csv")

### Eigenvalues ----
ev <- mfa_all$eig[,2]

saveRDS(ev, "data/output/mfa_all_ev.rds")

### Plot ----
#### Per elevation ----
pmfa_all_spp_ele <- ggplot(mfa_all_ind, aes(x = Dim.1, y = Dim.2, 
                                          color = species)) +
  geom_hline(yintercept = 0, color = "grey90", linewidth = 1) +     
  geom_vline(xintercept = 0, color = "grey90", linewidth = 1) + 
  geom_point(aes(shape = as.factor(elevation_m_asl)), size = 3, alpha = 0.9) +
  stat_ellipse(aes(group = species, colour = species), linewidth = 0.8) +
  geom_segment(data = mfa_all_traits,
               aes(x = 0, y = 0, xend = Dim1_scaled_12, yend = Dim2_scaled_12),
               inherit.aes = FALSE,
               arrow = arrow(length = unit(0.18, "cm")),
               color = "gray30", linewidth = 0.6) +
  geom_text_repel(
    data = mfa_all_traits, inherit.aes = FALSE,
    aes(x = Dim1_scaled_12, y = Dim2_scaled_12, label = traits),
    color = "gray20",  
    size = 5,
    max.overlaps = Inf,
    force = 1,           
    force_pull = 1,      
    min.segment.length = 0,  
    segment.size = 0.3,  
    box.padding = 0.5,
    show.legend = FALSE,
    fontface = "bold") +
  coord_equal() +
  scale_colour_manual(values = species_colors, name = "Species",
                      labels = species_labels) +
  scale_shape_manual(values = elevation_shapes, name = "Elevation (masl)") +
  labs(x = sprintf("Dim 1 (%.1f%%)", ev[1]),
       y = sprintf("Dim 2 (%.1f%%)", ev[2])) +
  theme_bw(base_size = 16) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.ticks = element_line(linewidth = 1),
    axis.line = element_line(linewidth = 1, colour = "black"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(1,1,1,1, "mm"),
    aspect.ratio = 1)

pmfa_all_spp_ele

##### Save plots ----
ggsave("results/img/pmfa_all_spp_ele.png", pmfa_all_spp_ele,
       width = 20, height = 15, units = "cm", dpi = 300)

#### Per species ----
pmfa_all_spp_grp <- ggplot(mfa_all_ind, aes(x = Dim.1, y = Dim.2)) +
  geom_hline(yintercept = 0, color = "grey90", linewidth = 1) +     
  geom_vline(xintercept = 0, color = "grey90", linewidth = 1) + 
  geom_point(aes(shape = species), 
             size = 3, color = "grey70") +
  stat_ellipse(aes(linetype = species), 
               linewidth = 0.8, show.legend = TRUE, color = "grey30") +
  geom_segment(data = mfa_all_traits,
               aes(x = 0, xend = Dim1_scaled_12,y = 0, yend = Dim2_scaled_12, color = group),
               inherit.aes = FALSE,
               arrow = arrow(length = unit(0.18, "cm")),
               linewidth = 1.2) +
  geom_text_repel(
    data = mfa_all_traits, inherit.aes = FALSE,
    aes(x = Dim1_scaled_12, y = Dim2_scaled_12, label = traits, color = group),
    size = 5,
    max.overlaps = Inf,
    force = 1,           
    force_pull = 1,      
    min.segment.length = 0,  
    segment.size = 0.3,  
    box.padding = 0.5,
    show.legend = FALSE,
    fontface = "bold") +
  coord_equal() +
  scale_shape_manual(values = species_shapes, name = "Species",
                     labels = species_labels) +
  scale_linetype_manual(values = species_ellipses, 
                        name = "Species",
                        labels = species_labels) +
  scale_color_manual(values = group_colors, name = "Trait group") +
  guides(
    shape = guide_legend(
      title = "Species",
      override.aes = list(
        linetype = species_ellipses,
        size = 4,
        color = "grey30"),
      ncol = 1,
      byrow = TRUE,
      keywidth = unit(1.5, "cm"),
      keyheight = unit(0.8, "cm")),
    linetype = "none",
    color = guide_legend(
      title = "Trait group",
      override.aes = list(
        linetype = "solid",
        alpha = 1,
        size = 1))) +
  labs(x = sprintf("Dim 1 (%.1f%%)", ev[1]),
       y = sprintf("Dim 2 (%.1f%%)", ev[2])) +
  theme_bw(base_size = 16) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.ticks = element_line(linewidth = 1),
    axis.line = element_line(linewidth = 1, colour = "black"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(1,1,1,1, "mm"),
    aspect.ratio = 1)

pmfa_all_spp_grp

##### Save plots ----
ggsave("results/img/pmfa_all_spp_grp.png", pmfa_all_spp_grp,
       width = 20, height = 15, units = "cm", dpi = 300)

#### Per Spp & Ele Dim 12 ----
pmfa_all_spp_ele_dim12 <- ggplot(mfa_all_ind, aes(x = Dim.1, y = Dim.2)) +
  geom_hline(yintercept = 0, color = "grey90", linewidth = 1) +     
  geom_vline(xintercept = 0, color = "grey90", linewidth = 1) + 
  geom_point(aes(shape = as.factor(elevation_m_asl)), 
             size = 3, color = "grey70") +
  stat_ellipse(aes(linetype = species), 
               linewidth = 0.8, show.legend = TRUE, color = "grey30") +
  geom_segment(data = mfa_all_traits,
               aes(x = 0, xend = Dim1_scaled_12,y = 0, yend = Dim2_scaled_12, color = group),
               inherit.aes = FALSE,
               arrow = arrow(length = unit(0.18, "cm")),
               linewidth = 1.2) +
  geom_text_repel(
    data = mfa_all_traits, inherit.aes = FALSE,
    aes(x = Dim1_scaled_12, y = Dim2_scaled_12, label = traits, color = group),
    size = 5,
    max.overlaps = Inf,
    force = 1,           
    force_pull = 1,      
    min.segment.length = 0,  
    segment.size = 0.3,  
    box.padding = 0.5,
    show.legend = FALSE,
    fontface = "bold") +
  coord_equal() +
  scale_shape_manual(values = elevation_shapes, name = "Elevation (m asl)") +
  scale_linetype_manual(values = species_ellipses, 
                        name = "Species",
                        labels = species_labels) +
  scale_color_manual(values = group_colors, name = "Trait group") +
  guides(
    shape = guide_legend(
      override.aes = list(linetype = "blank"),
      order = 1),
    linetype = guide_legend(
      override.aes = list(shape = NA),
      order = 2),
    color = guide_legend(order = 3)) +
  labs(x = sprintf("Dim 1 (%.1f%%)", ev[1]),
       y = sprintf("Dim 2 (%.1f%%)", ev[2])) +
  theme_bw(base_size = 16) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.ticks = element_line(linewidth = 1),
    axis.line = element_line(linewidth = 1, colour = "black"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(1,1,1,1, "mm"),
    aspect.ratio = 1)

pmfa_all_spp_ele_dim12

#### Per Spp & Ele Dim 13 ----
pmfa_all_spp_ele_dim13 <- ggplot(mfa_all_ind, aes(x = Dim.1, y = Dim.3)) +
  geom_hline(yintercept = 0, color = "grey90", linewidth = 1) +     
  geom_vline(xintercept = 0, color = "grey90", linewidth = 1) + 
  geom_point(aes(shape = as.factor(elevation_m_asl)), 
             size = 3, color = "grey70") +
  stat_ellipse(aes(linetype = species), 
               linewidth = 0.8, show.legend = TRUE, color = "grey30") +
  geom_segment(data = mfa_all_traits,
               aes(x = 0, xend = Dim1_scaled_13,y = 0, yend = Dim3_scaled_13, color = group),
               inherit.aes = FALSE,
               arrow = arrow(length = unit(0.18, "cm")),
               linewidth = 1.2) +
  geom_text_repel(
    data = mfa_all_traits, inherit.aes = FALSE,
    aes(x = Dim1_scaled_13, y = Dim3_scaled_13, label = traits, color = group),
    size = 5,
    max.overlaps = Inf,
    force = 1,           
    force_pull = 1,      
    min.segment.length = 0,  
    segment.size = 0.3,  
    box.padding = 0.5,
    show.legend = FALSE,
    fontface = "bold") +
  coord_equal() +
  scale_shape_manual(values = elevation_shapes, name = "Elevation (m asl)") +
  scale_linetype_manual(values = species_ellipses, 
                        name = "Species",
                        labels = species_labels) +
  scale_color_manual(values = group_colors, name = "Trait group") +
  guides(
    shape = guide_legend(
      override.aes = list(linetype = "blank"),
      order = 1),
    linetype = guide_legend(
      override.aes = list(shape = NA),
      order = 2),
    color = guide_legend(order = 3)) +
  labs(x = sprintf("Dim 1 (%.1f%%)", ev[1]),
       y = sprintf("Dim 3 (%.1f%%)", ev[3])) +
  theme_bw(base_size = 16) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.ticks = element_line(linewidth = 1),
    axis.line = element_line(linewidth = 1, colour = "black"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(1,1,1,1, "mm"),
    aspect.ratio = 1)

pmfa_all_spp_ele_dim13

#### Per Spp & Ele Dim 23 ----
pmfa_all_spp_ele_dim23 <- ggplot(mfa_all_ind, aes(x = Dim.2, y = Dim.3)) +
  geom_hline(yintercept = 0, color = "grey90", linewidth = 1) +     
  geom_vline(xintercept = 0, color = "grey90", linewidth = 1) + 
  geom_point(aes(shape = as.factor(elevation_m_asl)), 
             size = 3, color = "grey70") +
  stat_ellipse(aes(linetype = species), 
               linewidth = 0.8, show.legend = TRUE, color = "grey30") +
  geom_segment(data = mfa_all_traits,
               aes(x = 0, xend = Dim2_scaled_23,y = 0, yend = Dim3_scaled_23, color = group),
               inherit.aes = FALSE,
               arrow = arrow(length = unit(0.18, "cm")),
               linewidth = 1.2) +
  geom_text_repel(
    data = mfa_all_traits, inherit.aes = FALSE,
    aes(x = Dim2_scaled_23, y = Dim3_scaled_23, label = traits, color = group),
    size = 5,
    max.overlaps = Inf,
    force = 1,           
    force_pull = 1,      
    min.segment.length = 0,  
    segment.size = 0.3,  
    box.padding = 0.5,
    show.legend = FALSE,
    fontface = "bold") +
  coord_equal() +
  scale_shape_manual(values = elevation_shapes, name = "Elevation (m asl)") +
  scale_linetype_manual(values = species_ellipses, 
                        name = "Species",
                        labels = species_labels) +
  scale_color_manual(values = group_colors, name = "Trait group") +
  guides(
    shape = guide_legend(
      override.aes = list(linetype = "blank"),
      order = 1),
    linetype = guide_legend(
      override.aes = list(shape = NA),
      order = 2),
    color = guide_legend(order = 3)) +
  labs(x = sprintf("Dim 2 (%.1f%%)", ev[2]),
       y = sprintf("Dim 3 (%.1f%%)", ev[3])) +
  theme_bw(base_size = 16) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.ticks = element_line(linewidth = 1),
    axis.line = element_line(linewidth = 1, colour = "black"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(1,1,1,1, "mm"),
    aspect.ratio = 1)

pmfa_all_spp_ele_dim23

##### Save plots ----
ggsave("results/img/pmfa_all_spp_ele_dim12.png", pmfa_all_spp_ele_dim12,
       width = 20, height = 15, units = "cm", dpi = 300)

### Combined plot for Dim12, Dim13, Dim23 ----
# Get individual legends
elevation_legend <- get_legend(pmfa_all_spp_ele_dim12 + 
                                 guides(linetype = "none", color = "none") +
                                 theme(legend.key.height = unit(0.5, "cm")))

species_legend <- get_legend(pmfa_all_spp_ele_dim12 + 
                               guides(shape = "none", color = "none") +
                               theme(legend.key.height = unit(0.5, "cm")))

trait_legend <- get_legend(pmfa_all_spp_ele_dim12 + 
                             guides(shape = "none", linetype = "none"))


legend_combined <- plot_grid(elevation_legend, species_legend, trait_legend,
                             ncol = 2,
                             align = "v", axis = "t")

mfa_all_combined_plot1 <- plot_grid(pmfa_all_spp_ele_dim12 + theme(legend.position="none"), 
                                   pmfa_all_spp_ele_dim13 + theme(legend.position="none"), 
                                   pmfa_all_spp_ele_dim23 + theme(legend.position="none"),
                                   legend_combined,
                                   align = c("hv")) +
  theme(plot.background = element_rect(fill = "white", color = NA))

mfa_all_combined_plot1

##### Save plots ----
ggsave("results/img/mfa_all_combined_plot2.png", mfa_all_combined_plot1,
       width = 30, height = 30, units = "cm", dpi = 300)

## PERMANOVA  ----
### Matrix of distances ----
mfa_all_distance <- dist(mfa_table |> 
                           select(starts_with("Dim")))
### Species alone ----
set.seed(1)
adonis2(mfa_all_distance ~ species, data = mfa_table, permutations = 4999)

# Multivarite homogeneity
mfa_all_sp_bd <- betadisper(mfa_all_distance, mfa_table$species)
anova(mfa_all_sp_bd)
permutest(mfa_all_sp_bd, 999)
plot(mfa_all_sp_bd)

#### Pairwise comparisons among species ----
set.seed(1)
pw_spp_results <- pairwise.adonis2(mfa_all_distance ~ species, 
                 data = as.data.frame(mfa_table),
                 permutations = 4999, 
                 p.adjust.m = "BH")

pw_ss_table <- map_dfr(pw_spp_results, \(x) {
  #Check if x is a data.frame/matrix
  if(!is.data.frame(x) && !is.matrix(x)) {
    return(tibble())
  }
  
  if (nrow(x) < 3) {
    return(tibble())
  }
  
  #Extract Model row which contains the statistics
  model_row <- x[1,]
  total_row <- x[3,]
  
  #Create a tibble with the relevant statistics
  tibble(
    F = model_row$F,
    R2 = model_row$R2,
    p.value = model_row$`Pr(>F)`,
    SS = model_row$SumOfSqs,
    DF = paste0(model_row$Df, "/", total_row$Df),
    significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      p.value < 0.1 ~ ".",
      TRUE ~ "ns"))},
  .id = "Comparison")

print(pw_ss_table)

### Elevation alone ----
set.seed(1)
adonis2(mfa_all_distance ~ elevation_m_asl, data = mfa_table, permutations = 4999)

# Multivarite homogeneity
mfa_all_ele_bd <- betadisper(mfa_all_distance, mfa_table$elevation_m_asl)
anova(mfa_all_ele_bd)
permutest(mfa_all_ele_bd, 4999)
plot(mfa_all_ele_bd)

#### Pairwise comparisons among species ----
set.seed(1)
pw_ele_results <- pairwise.adonis2(mfa_all_distance ~ elevation_m_asl, 
                                   data = as.data.frame(mfa_table),
                                   permutations = 4999, 
                                   p.adjust.m = "BH")

pw_ele_table <- map_dfr(pw_ele_results, \(x) {
  #Check if x is a data.frame/matrix
  if(!is.data.frame(x) && !is.matrix(x)) {
    return(tibble())
  }
  
  if (nrow(x) < 3) {
    return(tibble())
  }
  
  #Extract Model row which contains the statistics
  model_row <- x[1,]
  total_row <- x[3,]
  
  #Create a tibble with the relevant statistics
  tibble(
    F = model_row$F,
    R2 = model_row$R2,
    p.value = model_row$`Pr(>F)`,
    SS = model_row$SumOfSqs,
    DF = paste0(model_row$Df, "/", total_row$Df),
    significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      p.value < 0.1 ~ ".",
      TRUE ~ "ns"))},
  .id = "Comparison")

print(pw_ele_table)


### Dimensions as variables ----
#### Dim.1 ---- 
mfa_all_dim1 <- feols(Dim.1 ~ species, mfa_all_ind)
summary(mfa_all_dim1)
car::Anova(mfa_all_dim1, type = 2)

# Perform Tukey HSD using emmeans
tukey_all_dim1 <- emmeans(mfa_all_dim1, pairwise ~ species, adjust = "tukey")

# View the results
tukey_all_dim1$contrasts
tukey_all_dim1$emmeans

# Create a nice table
tukey_t_all_dim1<- tukey_all_dim1$contrasts |> 
  as.data.frame() |> 
  mutate(
    estimate = round(estimate, 2),
    SE = round(SE, 2),
    t.ratio = round(t.ratio, 2),
    p.value = round(p.value, 4),
    significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      p.value < 0.1 ~ ".",
      TRUE ~ "ns"
    )
  )

print(tukey_t_all_dim1)

# Get compact letter display
cld_all_dim1 <- cld(tukey_all_dim1$emmeans, alpha = 0.05, Letters = letters)

print(cld_all_dim1)

# Calculate individual y-positions for each species (more precise)
species_max_dim1 <- mfa_all_ind |> 
  group_by(species) |> 
  summarise(
    max_value = max(Dim.1, na.rm = TRUE),
    q75 = quantile(Dim.1, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) |> 
  mutate(
    letter_y = pmax(max_value, q75) + 0.3  # Position above max or Q3, whichever is higher
  )

# Merge with cld_results
letters_all_dim1 <- cld_all_dim1 |> 
  as.data.frame() |> 
  mutate(
    species = factor(species, levels = species_order),
    .group = str_trim(.group)
  ) |> 
  left_join(species_max_dim1, by = "species")

##### Plot ----
mfa_all_dim1_p <- ggplot(mfa_all_ind, aes(x = species, y = Dim.1, color = species)) +
  geom_point(aes(y = Dim.1)) +
  geom_boxplot(aes(fill = species), alpha = 0.5) +
  geom_text(data = letters_all_dim1, 
            aes(x = species, y = letter_y, label = .group),
            color = "black", size = 6, fontface = "bold",
            inherit.aes = FALSE) +
  scale_color_manual(values = species_colors, labels = species_labels, name = "Species") +
  scale_fill_manual(values = species_colors, labels = species_labels, , name = "Species") +
  scale_x_discrete(labels = species_labels) +
  labs(
    x = "",
    y = "MFA Dimension 1"
  ) +
  theme_bw(base_size = 16) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.ticks = element_line(linewidth = 1),
    axis.line = element_line(linewidth = 1, colour = "black"),
    # legend.title = element_text(face = "bold"),
    # legend.text = element_text(size = 16),
    plot.margin = margin(1,1,1,1, "mm"),
    legend.position = "none",
    aspect.ratio = NULL)

mfa_all_dim1_p

#### Dim.2 ---- 
mfa_all_dim2 <- feols(Dim.2 ~ species, mfa_all_ind)
summary(mfa_all_dim2)
etable(mfa_all_dim2)

# Perform Tukey HSD using emmeans
tukey_all_dim2 <- emmeans(mfa_all_dim2, pairwise ~ species, adjust = "tukey")

# View the results
tukey_all_dim2$contrasts
tukey_all_dim2$emmeans

# Create a nice table
tukey_t_all_dim2<- tukey_all_dim2$contrasts |> 
  as.data.frame() |> 
  mutate(
    estimate = round(estimate, 2),
    SE = round(SE, 2),
    t.ratio = round(t.ratio, 2),
    p.value = round(p.value, 4),
    significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      p.value < 0.1 ~ ".",
      TRUE ~ "ns"
    )
  )

print(tukey_t_all_dim2)

# Get compact letter display
cld_all_dim2 <- cld(tukey_all_dim2$emmeans, alpha = 0.05, Letters = letters)

print(cld_all_dim2)

# Calculate individual y-positions for each species (more precise)
species_max_dim2 <- mfa_all_ind |> 
  group_by(species) |> 
  summarise(
    max_value = max(Dim.2, na.rm = TRUE),
    q75 = quantile(Dim.2, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) |> 
  mutate(
    letter_y = pmax(max_value, q75) + 0.3  # Position above max or Q3, whichever is higher
  )

# Merge with cld_results
letters_all_dim2 <- cld_all_dim2 |> 
  as.data.frame() |> 
  mutate(
    species = factor(species, levels = species_order),
    .group = str_trim(.group)
  ) |> 
  left_join(species_max_dim2, by = "species")

##### Plot ----
mfa_all_dim2_p <- ggplot(mfa_all_ind, aes(x = species, y = Dim.2, color = species)) +
  geom_point(aes(y = Dim.2)) +
  geom_boxplot(aes(fill = species), alpha = 0.5) +
  geom_text(data = letters_all_dim2, 
            aes(x = species, y = letter_y, label = .group),
            color = "black", size = 6, fontface = "bold",
            inherit.aes = FALSE) +
  scale_color_manual(values = species_colors, labels = species_labels, name = "Species") +
  scale_fill_manual(values = species_colors, labels = species_labels, , name = "Species") +
  scale_x_discrete(labels = species_labels) +
  labs(
    x = "",
    y = "MFA Dimension 2"
  ) +
  theme_bw(base_size = 16) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.ticks = element_line(linewidth = 1),
    axis.line = element_line(linewidth = 1, colour = "black"),
    # legend.title = element_text(face = "bold"),
    # legend.text = element_text(size = 16),
    plot.margin = margin(1,1,1,1, "mm"),
    legend.position = "none",
    aspect.ratio = NULL)

mfa_all_dim2_p

##### Save combined plot ----

mfa_dim1_2_plot <- mfa_all_dim1_p + mfa_all_dim2_p +
  plot_annotation(tag_levels = list(c("A.", "B."))) +
  plot_layout(ncol = 2) &
  theme(plot.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(2, 2, 2, 2, "mm"),
        plot.tag = element_text(size = 16, face = "bold"))

ggsave("results/img/mfa_dim1_2_plot.png", mfa_dim1_2_plot, 
       width = 10, height = 5, dpi = 300)


## FactoExtra plots ----
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


# 4. MFA for each species ----
## Define function to process one species ----
process_species <- function(sp, trait_data_wide) {
  cat("Processing species:", sp, "\n")
  
  # Get data for current species (without species column)
  sp_data <- trait_data_wide |> 
    filter(species == sp) |> 
    select(id, leaf_thickness, ldmc, sla, 
           bi, rd, rdmc, rtd, srl, 
           veg_height, root_depth, bgb_agb) |> 
    rename(LT = leaf_thickness, LDMC = ldmc, SLA = sla, 
           BI = bi, RD = rd, RDMC = rdmc, RTD = rtd, SRL = srl, 
           VHeight = veg_height, RDepth = root_depth, `BG:AG` = bgb_agb)
  
  # Remove missing data but keep track of which rows remain
  sp_clean <- na.omit(sp_data)
  
  if (nrow(sp_clean) < 3) {
    cat("  Skipping", sp, "- insufficient data\n")
    return(NULL)
  }
  
  # Store the row IDs before running MFA
  kept_rowids <- sp_clean$id
  
  # Remove rowid column for MFA
  mfa_data <- sp_clean |> select(-id)
  
  # Calculate ncp for this species
  grp_sizes <- c(3, 5, 3)
  ncp_max_sp <- min(nrow(mfa_data) - 1L, sum(grp_sizes) - length(grp_sizes))
  
  # Run MFA
  result <- MFA(
    base = mfa_data,
    group = grp_sizes,
    type = c("s", "s", "s"),
    name.group = c("Leaf", "Roots", "Plant size"),
    ncp = ncp_max_sp,
    graph = FALSE
  )
  
  # Store IDs in multiple accessible places
  rownames(result$ind$coord) <- kept_rowids  # Add as rownames to coordinates
  
  cat("  Success! N =", nrow(mfa_data), "\n")
  return(result)
}

## Run the loop with function ----
mfa_by_species <- list()

for (sp in species_order) {
  result <- process_species(sp, trait_data_wide)
  if (!is.null(result)) {
    mfa_by_species[[sp]] <- result
  }
  rm(result, sp)
}


## MFA Contribution table ----
### Contribution table ----
spp_contrib <- imap_dfr(mfa_by_species, \(x, y) {
  spp_contri_matrix <- x$quanti.var$contrib |> 
    as.data.frame() |> 
    rownames_to_column(var = "traits") |> 
    mutate(species = y) |> 
  left_join(x = _, 
            y = x$summary.quanti |> 
              select(group, variable),
            by = join_by(traits == variable)) |> 
    mutate(
      group = case_when(
        group == 1 ~ "Leaf",
        group == 2 ~ "Roots",
        group == 3 ~ "Plant size")) |> 
    relocate(group, .after = "traits") |> 
    relocate(species, .before = "traits")
})


### Function to create plots for a specific dimension ----
create_dim_plots <- function(data, dimension, dim_name) {
  
  # Create plots in the specified order
  map(species_order, \(sp) {
    df <- data |> filter(species == sp)
    species_idx <- match(sp, species_order)
    
    # Calculate tag based on dimension and species position - ADD PERIOD HERE
    tag_letter <- if(dim_name == "Dim-1") {
      paste0(LETTERS[species_idx], ".")  # A., B., C., D., E. for Dim-1
    } else {
      paste0(LETTERS[species_idx + length(species_order)], ".")  # F., G., H., I., J. for Dim-2
    }
    
    df |> 
      arrange(desc(.data[[dimension]])) |> 
      mutate(traits = factor(traits, levels = traits)) |> 
      mutate(group = fct_relevel(group, "Leaf", "Roots", "Plant size")) |> 
      ggplot(aes(x = traits, y = .data[[dimension]])) +
      geom_col(aes(fill = group)) +
      scale_fill_manual(values = group_colors) +
      scale_y_continuous(limits = if(dim_name == "Dim-1") c(0, 25) else c(0,50)) +
      geom_hline(yintercept = mean(data[[dimension]]), linetype = "dashed") +
      labs(y = paste0("Contributions to ", dim_name, " (%)"),
           x = "",
           title = if(dim_name == "Dim-1") paste0(tag_letter, " ", sp) else paste0(tag_letter, " "),
           fill = "") +
      theme_classic(base_size = 14) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 14, face = "bold", hjust = 0), 
        legend.position = "none",
        plot.margin = margin(1, 1, 1, 1, "mm"),
        panel.spacing = unit(0.5, "lines"),
        aspect.ratio = 1)
  }) |> 
    set_names(species_order)
}

### Plots for both dimensions ----
dim1_plots <- create_dim_plots(spp_contrib, "Dim.1", "Dim-1")
dim2_plots <- create_dim_plots(spp_contrib, "Dim.2", "Dim-2")

### Combine plots: 5 columns (species) x 2 rows (dimensions) ----
spp_contrib_com_plot <- wrap_plots(
  c(dim1_plots, dim2_plots),
  ncol = 5,
  nrow = 2,
  byrow = TRUE) +  # Fill by columns first (all Dim1, then all Dim2)
  plot_layout(guides = "collect",
              tag_level = "new") &  # Collect legends
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 14, face = "bold"),
    plot.tag = element_text(size = 14, face = "bold"))

print(spp_contrib_com_plot)

### Save the plot ----
ggsave("results/img/spp_contrib_dim1_dim2.png", spp_contrib_com_plot, 
       width = 20, height = 10, dpi = 300)


## MFA results ----

### Ind. coords table ----
mfa_spp_ind <- trait_data_wide |>
  select(id, species, family, growth_form, elevation_m_asl) |> 
  left_join(x = _,
            y = imap_dfr(mfa_by_species, \(x, y) {
              spp_mfa_matrix <- x$ind$coord |> 
                as.data.frame() |> 
                rownames_to_column(var = "id") |> 
                mutate(species = y)}),
            by = join_by(id == id, species == species)) |> 
  mutate(species = factor(species, levels = species_order))  


write_csv(mfa_spp_ind, "data/output/mfa_spp_ind.csv")

### Traits loadings (arrows) + scale ----

mfa_spp_traits <- imap_dfr(mfa_by_species, \(mfa_obj, sp_name) {
  
  # Get coordinates
  ind <- mfa_obj$ind$coord |> 
    as_tibble()
  var <- mfa_obj$quanti.var$coord |> 
    as_tibble(rownames = "traits")
  
  # Get dimension names
  dim_names <- names(ind)
  n_traits <- nrow(var)
  
  # Create scaled dims for common pairs
  scaled_dims <- tibble(.rows = n_traits)
  
  # Dim 1 vs Dim 2 (most common)
  if (length(dim_names) >= 2) {
    rng_ind1 <- range(ind[[dim_names[1]]])
    rng_ind2 <- range(ind[[dim_names[2]]])
    rng_var1 <- range(var[[dim_names[1]]])
    rng_var2 <- range(var[[dim_names[2]]])
    
    sf_12 <- 0.9 * min(diff(rng_ind1)/diff(rng_var1),
                       diff(rng_ind2)/diff(rng_var2))
    
    scaled_dims <- scaled_dims |> 
      bind_cols(
        tibble(
          Dim1_scaled_12 = var[[dim_names[1]]] * sf_12,
          Dim2_scaled_12 = var[[dim_names[2]]] * sf_12
        )
      )
  }
  
  # Dim 1 vs Dim 3 (if exists)
  if (length(dim_names) >= 3) {
    rng_ind1 <- range(ind[[dim_names[1]]])
    rng_ind3 <- range(ind[[dim_names[3]]])
    rng_var1 <- range(var[[dim_names[1]]])
    rng_var3 <- range(var[[dim_names[3]]])
    
    sf_13 <- 0.9 * min(diff(rng_ind1)/diff(rng_var1),
                       diff(rng_ind3)/diff(rng_var3))
    
    scaled_dims <- scaled_dims |> 
      bind_cols(
        tibble(
          Dim1_scaled_13 = var[[dim_names[1]]] * sf_13,
          Dim3_scaled_13 = var[[dim_names[3]]] * sf_13
        )
      )
  }
  
  # Dim 2 vs Dim 3 (if exists)
  if (length(dim_names) >= 3) {
    rng_ind2 <- range(ind[[dim_names[2]]])
    rng_ind3 <- range(ind[[dim_names[3]]])
    rng_var2 <- range(var[[dim_names[2]]])
    rng_var3 <- range(var[[dim_names[3]]])
    
    sf_23 <- 0.9 * min(diff(rng_ind2)/diff(rng_var2),
                       diff(rng_ind3)/diff(rng_var3))
    
    scaled_dims <- scaled_dims |> 
      bind_cols(
        tibble(
          Dim2_scaled_23 = var[[dim_names[2]]] * sf_23,
          Dim3_scaled_23 = var[[dim_names[3]]] * sf_23
        )
      )
  }
  
  # Combine everything
  bind_cols(
    var |> mutate(species = sp_name),
    scaled_dims
  ) |> 
    left_join(
      mfa_obj$summary.quanti |> select(group, variable),
      by = join_by(traits == variable)
    ) |> 
    mutate(
      species = factor(species, levels = species_order),
      group = case_when(
        group == "1" ~ "Leaf",
        group == "2" ~ "Roots",
        group == "3" ~ "Plant size"),
      group = fct_relevel(group, "Roots", after = 1)) |> 
    relocate(group, .after = "traits") |> 
    relocate(species, .before = "traits")
})

write_csv(mfa_spp_traits, "data/output/mfa_spp_traits.csv")

### Plot ----
#### Extract eigenvalues for each species ----
eigenvalues_by_species <- map_dfr(mfa_by_species, \(mfa_obj) {
  ev <- mfa_obj$eig[,2]  # Get percentage of variance explained
  tibble(
    Dim.1 = ev[1],
    Dim.2 = ev[2],
    Dim.3 = if(length(ev) >= 3) ev[3] else NA,
    Dim.4 = if(length(ev) >= 4) ev[4] else NA,
    Dim.5 = if(length(ev) >= 5) ev[5] else NA,
    Dim.6 = if(length(ev) >= 6) ev[6] else NA,
    Dim.7 = if(length(ev) >= 7) ev[7] else NA,
    Dim.8 = if(length(ev) >= 5) ev[8] else NA,
  )}, .id = "species") |> 
  mutate(species = factor(species, levels = species_order))

write_rds(eigenvalues_by_species, "data/output/mfa_spp_ev.rds")

species_labels_2 <- eigenvalues_by_species |> 
  mutate(
    label = sprintf("%s\nDim 1: %.1f%% | Dim 2: %.1f%%", 
                    species, Dim.1, Dim.2)) |> 
  select(species, label) |>  # Two columns for deframe()
  deframe()


mfa_spp_plot12 <- ggplot(mfa_spp_ind |> 
                           mutate(species = factor(species, levels = species_order)), aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(shape = as.factor(elevation_m_asl)), size = 3, alpha = 0.9) +
  # stat_ellipse(aes(group = as.factor(elevation_m_asl), colour = as.factor(elevation_m_asl)), size = 0.8) +
  geom_segment(data = mfa_spp_traits,
               aes(x = 0, y = 0, xend = Dim1_scaled_12, yend = Dim2_scaled_12,
                   color = group),
               inherit.aes = FALSE,
               arrow = arrow(length = unit(0.18, "cm")),
               linewidth = 0.6) +
  geom_text_repel(
    data = mfa_spp_traits, inherit.aes = FALSE,
    aes(x = Dim1_scaled_12, y = Dim2_scaled_12, label = traits, color = group),
    size = 7,
    max.overlaps = Inf,
    force = 1,           
    force_pull = 1,      
    min.segment.length = 0,  
    segment.size = 0.3,  
    box.padding = 0.5,
    show.legend = FALSE) +
  coord_equal() +
  facet_wrap(~ factor(species, levels = species_order), labeller = labeller(.cols = species_labels_2)) +
  scale_shape_manual(values = elevation_shapes, name = "Elevation (m asl)") +
  scale_color_manual(values = group_colors, name = "Trait group") +
  labs(x = "Dimension 1", y = "Dimension 2") +
  theme_bw(base_size = 16) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.ticks = element_line(linewidth = 1),
    axis.line = element_line(linewidth = 1, colour = "black"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 16),
    strip.text = element_text(face = "bold"),
    plot.margin = margin(1,1,1,1, "mm"),
    aspect.ratio = 1)


### Save the plot ----
ggsave("results/img/mfa_spp_plot12.png", mfa_spp_plot12, 
       width = 14, height = 10, dpi = 300)

## PERMANOVA  ----
### Matrix of distances ----

adonis_mfa_spp <- mfa_spp_ind |> 
  nest(.by = species) |> 
  mutate(dist_matrix = lapply(data, function(nest_data) 
    dist(nest_data |> 
           select(starts_with("Dim"))))) |> 
  mutate(
    #Run adonis model for each species to evaluate if there are differences among elevations
    adonis_models = map2(dist_matrix, data, \(x, y) {
      set.seed(12345)
      adonis2(x ~ elevation_m_asl, data = y, permutations = 4999)}),
    #Run betadisper to test for homogenity of variances among elevations
    betadisper = map2(dist_matrix, data, \(x, y) {
      bd <- betadisper(x, y$elevation_m_asl)
      set.seed(12345)
      perm.test <- permutest(bd, permutations = 999)
            test <- perm.test$tab}),
    #Run pairwise comparisons 
    pairwise = map2(dist_matrix, data, \(x, y) {
      set.seed(12345)
      pairwise.adonis2(x ~ elevation_m_asl, data = y,
                       permutations = 4999, 
                       p.adjust.m = "BH")}),
    #Extract adonis results
    adonis_summary = map(adonis_models, \(x) {
      table <- tibble(x)
      tibble(
        F = table$F[1],
        DF = paste0(table$Df[1],"/",table$Df[3]),
        R2 = table$R2[1],
        p = table$"Pr(>F)"[1])}),
    #Extract betadisper results
    betadisper_summary = map(betadisper, \(x) {
      table = tibble(x)
      tibble(
        F_bd = table$F[1],
        DF_bd = paste0(table$Df[1],"/", table$Df[2]),
        P_bd = table$"Pr(>F)"[1])}),
  #Extract pairwise results
  pairwise_results = map(pairwise, \(x) {
    as.data.frame(x) |> 
      rownames_to_column() |>  
      pivot_longer(
        cols = -c(parent_call,rowname),
        names_to = c("Comparison", ".value"),
        names_pattern = "([^\\.]+)\\.(.*)")})) |> 
  unnest(c(adonis_summary, betadisper_summary)) |> 
  relocate(pairwise_results,data, dist_matrix, adonis_models, betadisper, pairwise, .after = last_col())

adonis_mfa_spp_2 <- adonis_mfa_spp |> 
  select(1:8)
