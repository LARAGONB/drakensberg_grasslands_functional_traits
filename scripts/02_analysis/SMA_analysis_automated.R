run_sma_multitraits <- function(data, species_col, trait_pairs, 
                                nperm = 1000, plot_points = TRUE) {
  library(tidyverse)
  library(lmodel2)
  library(ggplot2)
  
  # Helper function to extract SMA information
  extract_sma <- function(model) {
    rr <- model$regression.results
    sma_row_idx <- which(grepl("Standard", rr[, "Method"], ignore.case = TRUE))
    if (length(sma_row_idx) == 0) sma_row_idx <- 3
    
    slope  <- as.numeric(rr[sma_row_idx, "Slope"])
    intercept <- as.numeric(rr[sma_row_idx, "Intercept"])
    
    # Significance from permutation p-value on r
    pval <- model$P.param
    sig <- ifelse(pval < 0.05, "significant", "nonsignificant")
    
    tibble(
      slope = slope,
      intercept = intercept,
      p_value = pval,
      significance = sig
    )
  }
  
  # Store results
  all_results <- list()
  all_plots <- list()
  all_stats <- tibble()
  
  # Loop through trait pairs
  for (pair in trait_pairs) {
    
    t1 <- pair[1]
    t2 <- pair[2]
    
    message("Running SMA for: ", t1, " ~ ", t2)
    
    df <- data %>%
      select(all_of(species_col), all_of(t1), all_of(t2)) %>%
      drop_na()
    
    colnames(df) <- c("species", "x", "y")
    
    # ----------------------
    # Global SMA
    # ----------------------
    mod_global <- lmodel2(y ~ x, data = df, nperm = nperm)
    global_sma <- extract_sma(mod_global) %>%
      mutate(species = "Global",
             trait_x = t1,
             trait_y = t2)
    
    # ----------------------
    # Species SMA
    # ----------------------
    species_sma <- df %>%
      group_by(species) %>%
      group_modify(~ {
        mod <- lmodel2(y ~ x, data = .x, nperm = nperm)
        extract_sma(mod)
      }) %>%
      ungroup() %>%
      mutate(trait_x = t1, trait_y = t2)
    
    # Store stats
    all_stats <- bind_rows(all_stats, global_sma, species_sma)
    
    # ----------------------
    # Build fitted lines
    # ----------------------
    x_seq <- seq(min(df$x), max(df$x), length.out = 200)
    
    global_line <- global_sma %>%
      tidyr::crossing(x = x_seq) %>%
      mutate(y = intercept + slope * x)
    
    species_lines <- species_sma %>%
      tidyr::crossing(x = x_seq) %>%
      mutate(y = intercept + slope * x)
    
    # ----------------------
    # Plot
    # ----------------------
    p <- ggplot(df, aes(x = x, y = y))
    
    if (plot_points) {
      p <- p + geom_point(aes(color = species), alpha = 0.55)
    }
    
    p <- p +
      geom_line(
        data = species_lines,
        aes(
          x = x,
          y = y,
          color = species,
          linetype = significance,
          linewidth = significance
        )
      ) +
      geom_line(
        data = global_line,
        aes(
          x = x,
          y = y,
          linetype = significance,
          linewidth = significance
        ),
        color = "black"
      ) +
      scale_linetype_manual(values = c(significant="solid", nonsignificant="dashed")) +
      scale_linewidth_manual(values = c(significant=1.2, nonsignificant=0.8)) +
      theme_minimal(base_size = 14) +
      labs(
        # title = paste("SMA for", t1, "vs", t2),
        # subtitle = "Bold = significant SMA; dashed = non-significant",
        x = t1,
        y = t2,
        color = "Species",
        linetype = "Significance",
        linewidth = "Significance"
      )
    
    # Store plot
    plot_name <- paste0(t1, "_vs_", t2)
    all_plots[[plot_name]] <- p
    
    # Store model results
    all_results[[plot_name]] <- list(
      global = mod_global,
      species = species_sma,
      plot = p
    )
  }
  
  # Return everything
  return(list(
    stats_table = all_stats,
    plots = all_plots,
    models = all_results
  ))
}

trait_data_wide <- read_csv("data/processed/v_PFCT7_clean_functional_traits_2023_wide.csv")

trait_list <- list(
  # c("srl", "sla"),
  # c("rd", "leaf_thickness"),
  # c("rtd", "ldmc"),
  c('sla', 'leaf_thickness'),
  c('sla', 'ldmc'),
  c('leaf_thickness', 'ldmc'),
  c('srl', 'rd'),
  c('srl', 'rtd'),
  c('rd', 'rtd')
)

results <- run_sma_multitraits(
  data = trait_data_wide,
  species_col = "species",
  trait_pairs = trait_list,
  nperm = 0
)
results$stats_table
results$models

library(patchwork)

combo_plot <- results$plots[[1]] + results$plots[[2]] + results$plots[[3]] + results$plots[[4]] + results$plots[[5]] + results$plots[[6]] +
  plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = 'a')

combo_plot

ggsave(combo_plot, filename = 'results/img/sma_trait_trait.png',
       width = 13.70, height = , 4.41, dpi = 320)
