################################################################################
# RV analyses
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

# 2. Load data and create groups ----
mfa_all_ind <- read_csv("data/output/mfa_all_ind.csv")
mfa_all_traits <- read_csv("data/output/mfa_all_traits.csv")
mfa_all_ev <- read_rds("data/output/mfa_all_ev.rds")

mfa_spp_ind <- read_csv("data/output/mfa_spp_ind.csv")
mfa_spp_traits <- read_csv("data/output/mfa_spp_traits.csv")
mfa_spp_ev <- read_rds("data/output/mfa_spp_ev.rds")


leaf_vars <- mfa_all_traits |> 
  select(traits, group) |> 
  filter(group %in% "Leaf") |> 
  pull(traits)

root_vars <- mfa_all_traits |> 
  select(traits, group) |> 
  filter(group %in% "Roots") |> 
  pull(traits)

plant_size_vars <- mfa_all_traits |> 
  select(traits, group) |> 
  filter(group %in% "Plant size") |> 
  pull(traits)

all_vars <- c(leaf_vars, root_vars, plant_size_vars)

# 3. Functions to analyze data ----
## Rv function of variable-coordinate matrices ----
rv_coords <- function(A, B) {
  stopifnot(is.matrix(A), is.matrix(B))
  
  # Get dimensions
  nA <- nrow(A)
  nB <- nrow(B)
  nmax <- max(nA, nB)
  
  # Pad matrices with zeros if needed (to make them same size)
  if (nA < nmax) {
    A <- rbind(A, matrix(0, nrow = nmax - nA, ncol = ncol(A)))
  }
  if (nB < nmax) {
    B <- rbind(B, matrix(0, nrow = nmax - nB, ncol = ncol(B)))
  }
  
  Sa <- A %*% t(A)
  Sb <- B %*% t(B)
  num <- sum(Sa * Sb)
  den <- sqrt(sum(Sa * Sa) * sum(Sb *Sb))
  num/den
}

## Permutation test on MFA varibal coordinate (single pair) ----
rv_perm_on_varcoords <- function(varcoords, all_vars, group1_vars, group2_vars,
                                 group1_name = "Group 1", group2_name = "Group 2",
                                 nperm = 10000, seed = NULL, verbose = TRUE) {
  if (!is.null(seed)) set.seed(seed)
  
  p <- length(group1_vars)
  q <- length(group2_vars)
  if (p + q > length(all_vars)) stop("group sizes larger than available variables")
  
  stopifnot(all(group1_vars %in% rownames(varcoords)),
            all(group2_vars %in% rownames(varcoords)),
            all(all_vars %in% rownames(varcoords)))
  
  # Ensure varcoords is a matrix
  varcoords <- as.matrix(varcoords)
  
  # Observed RV
  A_obs <- as.matrix(varcoords[group1_vars, , drop = FALSE])
  B_obs <- as.matrix(varcoords[group2_vars, , drop = FALSE])
  obs_rv <- rv_coords(A_obs, B_obs)
  
  # Permutation
  perm_rv <- numeric(nperm)
  for (i in seq_len(nperm)) {
    g1 <- sample(all_vars, size = p, replace = FALSE)
    remaining <- setdiff(all_vars, g1)
    g2 <- sample(remaining, size = q, replace = FALSE)
    Aperm <- varcoords[g1, , drop = FALSE]
    Bperm <- varcoords[g2, , drop = FALSE]
    perm_rv[i] <- rv_coords(Aperm, Bperm)
  }
  
  p_left_tail  <- (sum(perm_rv <= obs_rv) + 1) / (nperm + 1) # Orthogonality test
  p_right_tail <- (sum(perm_rv >= obs_rv) + 1) / (nperm + 1) # Coupling test
  
  list(obs_rv = obs_rv, perm_rv = perm_rv, 
       p_left_tail = p_left_tail,
       p_right_tail = p_right_tail,
       p = p, q = q, group1 = group1_vars, group2 = group2_vars,
       group1_name = group1_name, group2_name = group2_name)
}

## Plot single permutation result ----
plot_perm <- function(res, title = NULL) {
  dfnull <- data.frame(rv = res$perm_rv)
  
  # Find max y-value and round up to nearest 0.5
  max_density <- max(ggplot_build(
    ggplot(dfnull, aes(x = rv)) +
      geom_histogram(aes(y = after_stat(density)), bins = 60))$data[[1]]$y)
  max_y <- ceiling(max_density * 2) / 2  # Round up to nearest 0.5
  
  # Find max x-value and round up to nearest 0.5
  max_rv <- max(ggplot_build(
    ggplot(dfnull, aes(x = rv)) +
      geom_histogram(aes(y = after_stat(density)), bins = 60))$data[[1]]$x)
  max_x <- ceiling(max_rv * 2) / 2  # Round up to nearest 0.5
  
  # Calculate text x-position with offset
  text_x <- max_x * 0.85      # 85% across the plot
  text_y <- max_y * 0.95      # 95% up the plot
  
  # Plot
  ggplot(dfnull, aes(x = rv)) +
    geom_histogram(aes(y = after_stat(density)), bins = 60, fill = "grey50", color = "grey50") +
    geom_vline(xintercept = res$obs_rv, color = "red", linewidth = 1) +
    annotate("text", x = text_x, y = text_y,
             vjust = 1, hjust = 1,
             label = sprintf("obs = %.2f\nP = %.4f", res$obs_rv, res$p_left_tail),
             color = "red", size = 5) +
    scale_y_continuous(breaks = seq(0, max_y, by = 1)) +
    scale_x_continuous(breaks = seq(0, max_x, by = 0.1)) +
    labs(x = paste("Absolute value of correlation\n", 
                   res$group1_name, 'and', res$group2_name),
         y = "Density", title = title) +
    theme_minimal(base_size = 16) +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 16, color = "black"),
      axis.ticks = element_line(linewidth = 1),
      axis.line = element_line(linewidth = 1, colour = "black"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 16),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(1,1,1,1, "mm"))
}


# MFA ALL ----
# 4. Analyses for MFA across all individual ----
varcoords <- mfa_all_traits |> 
  select(traits, starts_with("Dim.")) |> 
  column_to_rownames(var = "traits")

results_all_species <- list(
    leaf_vs_root = rv_perm_on_varcoords(varcoords, all_vars, leaf_vars, root_vars, "Leaf", "Roots", nperm = 10000, seed = 42),
    leaf_vs_plant_size = rv_perm_on_varcoords(varcoords, all_vars, leaf_vars, plant_size_vars, "Leaf", "Plant size", nperm = 10000, seed = 42),
    root_vs_plant_size = rv_perm_on_varcoords(varcoords, all_vars, root_vars, plant_size_vars, "Roots", "Plant size", nperm = 10000, seed = 42))

summary_all_species <- tibble(
  Species = c("All", "All", "All"),
  Comparison = c("Leaf_Roots", "Leaf_Plant_size", "Roots_Plant_size"),
  RV = c(results_all_species$leaf_vs_root$obs_rv, 
         results_all_species$leaf_vs_plant_size$obs_rv, 
         results_all_species$root_vs_plant_size$obs_rv),
  P_orthogonal = c(results_all_species$leaf_vs_root$p_left_tail, 
                   results_all_species$leaf_vs_plant_size$p_left_tail, 
                   results_all_species$root_vs_plant_size$p_left_tail),
  P_coupled = c(results_all_species$leaf_vs_root$p_right_tail, 
                results_all_species$leaf_vs_plant_size$p_right_tail, 
                results_all_species$root_vs_plant_size$p_right_tail))

print(summary_all_species)


# 5. Plots ----
plots_all_species <- lapply(results_all_species, function(p) {
  plot_perm(p)})

## A. Leaf vs roots ----
pres_leaf_roots <- plots_all_species$leaf_vs_root + 
  labs(y = "", x = "", title = "B. Leaf and Roots") +
  theme(
    plot.title = element_text(size = 14, face = "bold"))

## B. Leaf vs Plant size ----
pres_leaf_plant_size <- plots_all_species$leaf_vs_plant_size + 
  labs(y = "Density", x = "", title = "C. Leaf and Plant size") +
  theme(
    plot.title = element_text(size = 14, face = "bold"))

## C. Roots vs Plant size ----
pres_roots_plant_size <- plots_all_species$root_vs_plant_size + 
  labs(y = "", x = "Absolute value of correlation", title = "D. Roots and Plant size") +
  theme(
    plot.title = element_text(size = 14, face = "bold"))


# 6. Combine RV plots ----
rv_all_plots <- plot_grid(pres_leaf_roots,
                      pres_leaf_plant_size,
                      pres_roots_plant_size,
                      ncol = 1)
rv_all_plots

# 7. Create full plot ----
pmfa_all_spp_ele_dim12 <- ggplot(mfa_all_ind, aes(x = Dim.1, y = Dim.2)) +
  geom_hline(yintercept = 0, color = "grey90", linewidth = 1) +     
  geom_vline(xintercept = 0, color = "grey90", linewidth = 1) + 
  geom_point(aes(shape = as.factor(elevation_m_asl)), 
             size = 3, color = "grey70") +
  stat_ellipse(aes(linetype = species), 
               linewidth = 0.8, show.legend = TRUE, color = "grey30", ) +
  geom_segment(data = mfa_all_traits,
               aes(x = 0, xend = Dim1_scaled_12,y = 0, yend = Dim2_scaled_12, color = group),
               inherit.aes = FALSE,
               arrow = arrow(length = unit(0.18, "cm")),
               linewidth = 1.2) +
  geom_text_repel(
    data = mfa_all_traits, inherit.aes = FALSE,
    aes(x = Dim1_scaled_12, y = Dim2_scaled_12, label = traits, color = group),
    size = 9,
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
  labs(x = sprintf("Dim 1 (%.1f%%)", mfa_all_ev[1]),
       y = sprintf("Dim 2 (%.1f%%)", mfa_all_ev[2])) +
  theme_bw(base_size = 24) +
  theme(
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 24, color = "black"),
    axis.ticks = element_line(linewidth = 1),
    axis.line = element_line(linewidth = 1, colour = "black"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(1,1,1,1, "mm"),
    aspect.ratio = 1)

pmfa_all_spp_ele_dim12

full_plot <- plot_grid(pmfa_all_spp_ele_dim12,
                       rv_all_plots,
                       ncol = 2,
                       rel_widths = c(4, 2.5),
                       nrow = 1,
                       rel_heights = c(1, 1),
                       labels = "A.",
                       label_size = 14,
                       align = "v",
                       axis = "l") +
  theme(plot.background = element_rect(fill = "white", color = NA))
full_plot

### Save the plot ----
ggsave("results/img/full_plot.png", full_plot, 
       width = 18, height = 8, dpi = 300)

# MFA BY SPP ----

# 8. Apply pairwise RV tests to each species MFA ----
species_names <- unique(mfa_spp_ind$species)

results_by_species <- lapply(species_names, function(sp) {
  res_mfa_sp <- mfa_spp_traits |> 
    filter(species == sp)
  varcoords_sp <- res_mfa_sp |> 
    select(traits, starts_with("Dim.")) |> 
    column_to_rownames(var = "traits")
  
  # Check variables exist in this species' MFA
  if (!all(c(leaf_vars, root_vars, plant_size_vars) %in% rownames(varcoords_sp))) {
    message("Species ", sp, ": missing some variables, skipping")
    return(NULL)
  }
  
  list(
    species = sp,
    leaf_vs_root = rv_perm_on_varcoords(varcoords_sp, all_vars, leaf_vars, root_vars, "Leaf", "Roots", nperm = 10000, seed = 42),
    leaf_vs_plant_size = rv_perm_on_varcoords(varcoords_sp, all_vars, leaf_vars, plant_size_vars, "Leaf", "Plant size", nperm = 10000, seed = 42),
    root_vs_plant_size = rv_perm_on_varcoords(varcoords_sp, all_vars, root_vars, plant_size_vars, "Roots", "Plant size", nperm = 10000, seed = 42)
  )
})

names(results_by_species) <- species_names

# Extract summary
summary_by_species <- results_by_species |> 
  discard(is.null) |> 
  map(\(x) { tibble(
    Species = x$species,
    Comparison = c("Leaf_Roots", "Leaf_Plant_size", "Roots_Plant_size"),
    RV = c(x$leaf_vs_root$obs_rv, x$leaf_vs_plant_size$obs_rv, x$root_vs_plant_size$obs_rv),
    P_orthogonal = c(x$leaf_vs_root$p_left_tail, x$leaf_vs_plant_size$p_left_tail, x$root_vs_plant_size$p_left_tail),
    P_coupled = c(x$leaf_vs_root$p_right_tail, x$leaf_vs_plant_size$p_right_tail, x$root_vs_plant_size$p_right_tail)
  )}) |> 
  list_rbind()

print(summary_by_species)

## Plots ----
plots_by_species <- lapply(results_by_species, function(p) {
  list(
    pres_spp_leaf_roots = plot_perm(p$leaf_vs_root),
    pres_spp_leaf_plant_size = plot_perm(p$leaf_vs_plant_size),
    pres_spp_roots_plant_size = plot_perm(p$root_vs_plant_size))
  })


combined_plots_by_species <- lapply(names(plots_by_species), function(sp_idx) {
  sp_plots <- plots_by_species[[sp_idx]]
  
  combined <- plot_grid(
    sp_plots$pres_spp_leaf_roots +
      labs(x = "", y = "", title = "B. Leaf and Roots") +
      theme(
        plot.title = element_text(size = 14, face = "bold")),
    sp_plots$pres_spp_leaf_plant_size +
      labs(x = "", title = "C. Leaf and Plant size") +
      theme(
        plot.title = element_text(size = 14, face = "bold")),
    sp_plots$pres_spp_roots_plant_size +
      labs(x = "Absolute value of correlation", y = "", title = "D. Roots and Plant size") +
      theme(
        plot.title = element_text(size = 14, face = "bold")),
    ncol = 1)
  
  return(combined)
})

names(combined_plots_by_species) <- species_names


#MFA plot by species 
pmfa_by_species <- lapply(species_names, function(sp) {
  mfa_spp_ind_sp <- mfa_spp_ind |> 
    filter(species == sp)
  mfa_spp_traits_sp <- mfa_spp_traits |> 
    filter(species == sp)
  mfa_spp_ev_sp <- mfa_spp_ev |> 
    filter(species == sp)
  
  # Get unique elevations for this species
  elevations_sp <- unique(mfa_spp_ind_sp$elevation_m_asl)
  
  # Filter shape/linetype vectors to match species elevations
  shapes_sp <- elevation_shapes[names(elevation_shapes) %in% as.character(elevations_sp)]
  ellipses_sp <- elevation_ellipses[names(elevation_ellipses) %in% as.character(elevations_sp)]
  labels_sp <- elevation_labels[names(elevation_labels) %in% as.character(elevations_sp)]
  
  
  ggplot(mfa_spp_ind_sp, aes(x = Dim.1, y = Dim.2)) +
    geom_hline(yintercept = 0, color = "grey90", linewidth = 1) +     
    geom_vline(xintercept = 0, color = "grey90", linewidth = 1) + 
    geom_point(aes(shape = as.factor(elevation_m_asl)), 
               size = 3, color = "grey70") +
    stat_ellipse(aes(linetype = as.factor(elevation_m_asl)), 
                 linewidth = 0.8, show.legend = TRUE, color = "grey30") +
    geom_segment(data = mfa_spp_traits_sp,
                 aes(x = 0, xend = Dim1_scaled_12,y = 0, yend = Dim2_scaled_12, color = group),
                 inherit.aes = FALSE,
                 arrow = arrow(length = unit(0.18, "cm")),
                 linewidth = 1.2) +
    geom_text_repel(
      data = mfa_spp_traits_sp, inherit.aes = FALSE,
      aes(x = Dim1_scaled_12, y = Dim2_scaled_12, label = traits, color = group),
      size = 9,
      max.overlaps = Inf,
      force = 1,           
      force_pull = 1,      
      min.segment.length = 0,  
      segment.size = 0.3,  
      box.padding = 0.5,
      show.legend = FALSE,
      fontface = "bold") +
    coord_equal() +
    scale_shape_manual(values = shapes_sp, name = "Elevation (m asl)",
                       labels = labels_sp) +
    scale_linetype_manual(values = ellipses_sp, 
                          name = "Elevation (m asl)",
                          labels = labels_sp) +
    scale_color_manual(values = group_colors, name = "Trait group") +
    guides(
      shape = guide_legend(
        title = "Elevation (m asl)",
        override.aes = list(
          linetype = as.character(ellipses_sp),
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
    labs(x = sprintf("Dim 1 (%.1f%%)", mfa_spp_ev_sp$Dim.1[[1]]),
         y = sprintf("Dim 2 (%.1f%%)", mfa_spp_ev_sp$Dim.2[[1]])) +
    theme_bw(base_size = 24) +
    theme(
      axis.title = element_text(size = 24),
      axis.text = element_text(size = 24, color = "black"),
      axis.ticks = element_line(linewidth = 1),
      axis.line = element_line(linewidth = 1, colour = "black"),
      legend.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 16),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(1,1,1,1, "mm"),
      aspect.ratio = 1)
  
})
names(pmfa_by_species) <- species_names

### ERCA Plot ----
full_plot_ERCA <- plot_grid(pmfa_by_species$`Eragrostis capensis`,
                       combined_plots_by_species$`Eragrostis capensis`,
                       ncol = 2,
                       rel_widths = c(4, 2.5),
                       nrow = 1,
                       rel_heights = c(1, 1),
                       labels = "A.",
                       label_size = 14,
                       align = "v",
                       axis = "l") +
  theme(plot.background = element_rect(fill = "white", color = NA))
full_plot_ERCA

#### Save the plot ----
ggsave("results/img/full_plot_ERCA.png", full_plot_ERCA, 
       width = 18, height = 8, dpi = 300)

### HAFA ----
full_plot_HAFA <- plot_grid(pmfa_by_species$`Harpochloa falx`,
                            combined_plots_by_species$`Harpochloa falx`,
                            ncol = 2,
                            rel_widths = c(4, 2.5),
                            nrow = 1,
                            rel_heights = c(1, 1),
                            labels = "A.",
                            label_size = 14,
                            align = "v",
                            axis = "l") +
  theme(plot.background = element_rect(fill = "white", color = NA))
full_plot_HAFA

#### Save the plot ----
ggsave("results/img/full_plot_HAFA.png", full_plot_HAFA, 
       width = 18, height = 8, dpi = 300)

### THTR ----
full_plot_THTR <- plot_grid(pmfa_by_species$`Themeda triandra`,
                            combined_plots_by_species$`Themeda triandra`,
                            ncol = 2,
                            rel_widths = c(4, 2.5),
                            nrow = 1,
                            rel_heights = c(1, 1),
                            labels = "A.",
                            label_size = 14,
                            align = "v",
                            axis = "l") +
  theme(plot.background = element_rect(fill = "white", color = NA))
full_plot_THTR

#### Save the plot ----
ggsave("results/img/full_plot_THTR.png", full_plot_THTR, 
       width = 18, height = 8, dpi = 300)


### HEPI ----
full_plot_HEPI <- plot_grid(pmfa_by_species$`Helichrysum pilosellum`,
                            combined_plots_by_species$`Helichrysum pilosellum`,
                            ncol = 2,
                            rel_widths = c(4, 2.5),
                            nrow = 1,
                            rel_heights = c(1, 1),
                            labels = "A.",
                            label_size = 14,
                            align = "v",
                            axis = "l") +
  theme(plot.background = element_rect(fill = "white", color = NA))
full_plot_HEPI

#### Save the plot ----
ggsave("results/img/full_plot_HEPI.png", full_plot_HEPI, 
       width = 18, height = 8, dpi = 300)


### SEGL ----
full_plot_SEGL <- plot_grid(pmfa_by_species$`Senecio glaberrimus`,
                            combined_plots_by_species$`Senecio glaberrimus`,
                            ncol = 2,
                            rel_widths = c(4, 2.5),
                            nrow = 1,
                            rel_heights = c(1, 1),
                            labels = "A.",
                            label_size = 14,
                            align = "v",
                            axis = "l") +
  theme(plot.background = element_rect(fill = "white", color = NA))
full_plot_SEGL

#### Save the plot ----
ggsave("results/img/full_plot_SEGL.png", full_plot_SEGL, 
       width = 18, height = 8, dpi = 300)


# 9.Full table summary ----

summary_rv_complete <- rbind(summary_all_species, summary_by_species)
summary_rv_complete |> 
  filter(P_orthogonal <= 0.1)
