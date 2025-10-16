# Rv function of variable-coordinate matrices
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

# Permutation test on MFA varibal coordinate (single pair)
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

# Plot single permutation result
plot_perm <- function(res, title = NULL) {
  dfnull <- data.frame(rv = res$perm_rv)
  
  # Calculate text x-position with offset
  x_range <- res$obs_rv - min(dfnull)
  text_x <- res$obs_rv - x_range
  
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
  
  ggplot(dfnull, aes(x = rv)) +
    geom_histogram(aes(y = after_stat(density)), bins = 60, fill = "grey50", color = "grey50") +
    geom_vline(xintercept = res$obs_rv, color = "red", linewidth = 1) +
    annotate("text", x = text_x, y = Inf,
             label = sprintf("obs = %.2f\nP = %.4f", res$obs_rv, res$p_left_tail),
             vjust = 2, hjust = 0,
             color = "red", size = 5) +
    scale_y_continuous(breaks = seq(0, max_y, by = 1)) +
    scale_x_continuous(breaks = seq(0, max_x, by = 0.1)) +
    labs(x = paste("Absolute value of correlation between", res$group1_name, 'and', res$group2_name),
         y = "Density", title = title) +
    theme_bw(base_size = 16) +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 16, color = "black"),
      axis.ticks = element_line(linewidth = 1),
      axis.line = element_line(linewidth = 1, colour = "black"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 16),
      # panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(1,1,1,1, "mm"))
}

varcoords <- mfa_all_traits |> 
  select(traits, starts_with("Dim.")) |> 
  column_to_rownames(var = "traits")


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

# 1. Leaf vs roots
res_leaf_roots <- rv_perm_on_varcoords(varcoords = varcoords, all_vars = all_vars,
                                       group1_vars = leaf_vars, group2_vars = root_vars,
                                       group1_name = "Leaf", group2_name = "Roots",
                                       nperm = 10000, seed = 42)
cat("Leaf vs Roots:\n")
cat("  Observed RV =", res_leaf_roots$obs_rv, "\n")
cat("  p (left-tail, orthogonality) =", res_leaf_roots$p_left_tail, "\n")
cat("  p (right-tail, coupling) =", res_leaf_roots$p_right_tail, "\n")

pres_leaf_roots <- plot_perm(res_leaf_roots)
pres_leaf_roots

# 1. Leaf vs plant_size
res_leaf_plant_size <- rv_perm_on_varcoords(varcoords = varcoords, all_vars = all_vars,
                                       group1_vars = leaf_vars, group2_vars = plant_size_vars,
                                       group1_name = "Leaf", group2_name = "Plant size",
                                       nperm = 10000, seed = 42)
cat("Leaf vs plant_size:\n")
cat("  Observed RV =", res_leaf_plant_size$obs_rv, "\n")
cat("  p (left-tail, orthogonality) =", res_leaf_plant_size$p_left_tail, "\n")
cat("  p (right-tail, coupling) =", res_leaf_plant_size$p_right_tail, "\n")

pres_leaf_plant_size <- plot_perm(res_leaf_plant_size)
pres_leaf_plant_size

# 1. roots vs plant_size
res_roots_plant_size <- rv_perm_on_varcoords(varcoords = varcoords, all_vars = all_vars,
                                            group1_vars = root_vars, group2_vars = plant_size_vars,
                                            group1_name = "Roots", group2_name = "Plant size",
                                            nperm = 10000, seed = 42)
cat("roots vs plant_size:\n")
cat("  Observed RV =", res_roots_plant_size$obs_rv, "\n")
cat("  p (left-tail, orthogonality) =", res_roots_plant_size$p_left_tail, "\n")
cat("  p (right-tail, coupling) =", res_roots_plant_size$p_right_tail, "\n")

pres_roots_plant_size <- plot_perm(res_roots_plant_size)
pres_roots_plant_size
