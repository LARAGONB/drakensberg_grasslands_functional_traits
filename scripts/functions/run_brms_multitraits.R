library(brms)
library(tidybayes)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

species_names <- c(
  "Eragrostis capensis",
  "Harpochloa falx",
  "Themeda triandra",
  "Helichrysum pilosellum",
  "Senecio glaberrimus")

species_labels <- c(
  "Eragrostis capensis" = "ERCA",
  "Harpochloa falx" = "HAFA",
  "Themeda triandra" = "THTR",
  "Helichrysum pilosellum" = "HEPI",
  "Senecio glaberrimus" = "SEGL")

# Main function to run models over multiple trait pairs
run_trait_models <- function(data, trait_list, group = "species", random_slope = TRUE,
                             chains = 4, cores = 4, ndraws = 1000) {
  
  # Function for a single trait pair
  fit_single_pair <- function(traits) {
    response <- traits[2]
    predictor <- traits[1]
    
    # Formula
    f <- if(random_slope){
      as.formula(paste0(response, " ~ ", predictor, " + (1 + ", predictor, " | ", group, ")"))
    } else {
      as.formula(paste0(response, " ~ ", predictor, " + (1 | ", group, ")"))
    }
    
    # Fit model
    fit <- brm(f, data = data, chains = chains, cores = cores)
    
    # Key brms summary table
    key_summary <- summary(fit)$fixed %>%
      as.data.frame() %>%
      mutate(parameter = rownames(summary(fit)$fixed),
             response = response,
             predictor = predictor)
    
    # Posterior samples
    post <- posterior_samples(fit)
    
    # Population-level slope
    pop_slope <- post %>%
      select(starts_with(paste0("b_", predictor))) %>%
      rename(slope = 1) %>%
      summarise(
        mean = mean(slope),
        lower = quantile(slope, 0.025),
        upper = quantile(slope, 0.975),
        prob_positive = mean(slope > 0),
        prob_negative = mean(slope < 0)
      ) %>%
      mutate(level = "Population", parameter = predictor, response = response)
    
    pop_post <- posterior_samples(fit, pars = paste0("b_", predictor))[[1]]
    pop_credible <- ifelse(quantile(pop_post, 0.025) <= 0 & quantile(pop_post, 0.975) >= 0,
                           "dashed", "solid")
    
    
    # Species-level slopes (if random slope)
    species_slope <- NULL
    if(random_slope){
      species_pars <- grep(paste0("r_", group, ".*", predictor), colnames(post), value = TRUE)
      species_slope <- post %>%
        select(all_of(c(paste0("b_", predictor), species_pars))) %>%
        mutate(draw = row_number()) %>%
        pivot_longer(-draw, names_to = "species_par", values_to = "value") %>%
        mutate(
          species = gsub(paste0("r_", group, "\\[(.*),", predictor, "\\]"), "\\1", species_par),
          slope = value
        ) %>%
        group_by(species) %>%
        summarise(
          mean = mean(slope),
          lower = quantile(slope, 0.025),
          upper = quantile(slope, 0.975),
          prob_positive = mean(slope > 0),
          prob_negative = mean(slope < 0),
          .groups = "drop"
        ) %>%
        mutate(level = "Species", parameter = predictor, response = response)
    }
    
    slope_table <- bind_rows(pop_slope, species_slope)
    
    # 1) Population-level prediction grid (over full predictor range)
    pred_global_grid <- tibble(
      !!predictor := seq(min(data[[predictor]], na.rm = TRUE),
                         max(data[[predictor]], na.rm = TRUE), length.out = 100)
    )
    
    pred_global <- add_epred_draws(
      pred_global_grid,
      object = fit,
      re_formula = NA
    )
    
    pred_global_summary <- pred_global %>%
      group_by(!!sym(predictor)) %>%
      mean_qi(.epred)
    pred_global_summary$linetype <- pop_credible
    
    # 2) Species-specific prediction grid (restricted to observed range)
    species_pars <- grep(paste0("r_", group, ".*", predictor), colnames(post), value = TRUE)
    
    species_credible <- sapply(species_pars, function(sp) {
      sp_post <- post[[paste0("b_", predictor)]] + post[[sp]]  # total slope
      if(quantile(sp_post, 0.025) <= 0 & quantile(sp_post, 0.975) >= 0) "dashed" else "solid"
    })
    
    species_linetypes <- data.frame(species = species_names, linetype = species_credible[2:6])
    
    species_ranges <- data %>%
      group_by(.data[[group]]) %>%
      summarise(
        min_pred = min(.data[[predictor]]),
        max_pred = max(.data[[predictor]]),
        .groups = "drop"
      )
    
    pred_species_grid <- species_ranges %>%
      rowwise() %>%
      mutate(
        !!predictor := list(seq(min_pred, max_pred, length.out = 50))
      ) %>%
      unnest(cols = !!predictor)
    
    pred_species <- add_epred_draws(
      pred_species_grid,
      object = fit,
      re_formula = NULL
    )
    
    pred_species_summary <- pred_species %>%
      group_by(.data[[group]], !!sym(predictor)) %>%
      mean_qi(.epred) %>%
      left_join(species_linetypes, by = c("species"))
    
    # Plot with global line + CI and species lines
    p <- ggplot(data, aes_string(x = predictor, y = response)) +
      geom_point(alpha = 0.55, aes_string(color = group)) +
      
      # # Population-level 95% CI ribbon
      # stat_lineribbon(data = pred_global_summary,
      #             aes_string(x = predictor, ymin = ".lower", ymax = ".upper"),
      #             fill = "grey70", alpha = 0.4) +
      # Population-level line
      geom_line(data = pred_global_summary,
                aes_string(x = predictor, y = ".epred", linetype = 'linetype'),
                linewidth = 1.2, color = "black") +
      # Species lines
      geom_line(data = pred_species_summary,
                aes_string(x = predictor, y = ".epred", color = group, linetype = "linetype"),
                linewidth = 1) +

      scale_color_manual(values = c("#42049EFF","#8204A7FF","#B6308BFF","#F79143FF","#FCCE25FF"), labels = species_labels) +
      scale_linetype_identity() +
      theme_classic()
    
    p
    
    return(list(
      fit = fit,
      key_summary = key_summary,
      slope_table = slope_table,
      plot = p
    ))
  }
  
  # Loop over all trait pairs
  results <- map(trait_list, fit_single_pair)
  names(results) <- sapply(trait_list, function(x) paste(x[2], "vs", x[1], sep = "_"))
  
  return(results)
}


# Main function to run models over multiple trait pairs
run_trait_models2 <- function(data, trait_list, group = "species", random_slope = TRUE,
                             chains = 4, cores = 4, ndraws = 1000) {
  
  # Function for a single trait pair
  fit_single_pair <- function(traits) {
    response <- traits[2]
    predictor <- traits[1]
    
    # Formula
    f <- if(random_slope){
      as.formula(paste0(response, " ~ ", predictor, " + (1 + ", predictor, " | ", group, ")"))
    } else {
      as.formula(paste0(response, " ~ ", predictor, " + (1 | ", group, ")"))
    }
    
    # Fit model
    fit <- brm(f, data = data, chains = chains, cores = cores)
    
    # Key brms summary table
    key_summary <- summary(fit)$fixed %>%
      as.data.frame() %>%
      mutate(parameter = rownames(summary(fit)$fixed),
             response = response,
             predictor = predictor)
    
    # Posterior samples
    post <- posterior_samples(fit)
    
    # Population-level slope
    pop_slope <- post %>%
      select(starts_with(paste0("b_", predictor))) %>%
      rename(slope = 1) %>%
      summarise(
        mean = mean(slope),
        lower = quantile(slope, 0.025),
        upper = quantile(slope, 0.975),
        prob_positive = mean(slope > 0),
        prob_negative = mean(slope < 0)
      ) %>%
      mutate(level = "Population", parameter = predictor, response = response)
    
    pop_post <- posterior_samples(fit, pars = paste0("b_", predictor))[[1]]
    pop_credible <- ifelse(quantile(pop_post, 0.025) <= 0 & quantile(pop_post, 0.975) >= 0,
                           "dashed", "solid")
    
    
    # Species-level slopes (if random slope)
    species_slope <- NULL
    if(random_slope){
      species_pars <- grep(paste0("r_", group, ".*", predictor), colnames(post), value = TRUE)
      species_slope <- post %>%
        select(all_of(c(paste0("b_", predictor), species_pars))) %>%
        mutate(draw = row_number()) %>%
        pivot_longer(-draw, names_to = "species_par", values_to = "value") %>%
        mutate(
          species = gsub(paste0("r_", group, "\\[(.*),", predictor, "\\]"), "\\1", species_par),
          slope = value
        ) %>%
        group_by(species) %>%
        summarise(
          mean = mean(slope),
          lower = quantile(slope, 0.025),
          upper = quantile(slope, 0.975),
          prob_positive = mean(slope > 0),
          prob_negative = mean(slope < 0),
          .groups = "drop"
        ) %>%
        mutate(level = "Species", parameter = predictor, response = response)
    }
    
    slope_table <- bind_rows(pop_slope, species_slope)
    
    # 1) Population-level prediction grid (over full predictor range)
    pred_global_grid <- tibble(
      !!predictor := seq(min(data[[predictor]], na.rm = TRUE),
                         max(data[[predictor]], na.rm = TRUE), length.out = 100)
    )
    
    pred_global <- add_epred_draws(
      pred_global_grid,
      object = fit,
      re_formula = NA
    )
    
    pred_global_summary <- pred_global %>%
      group_by(!!sym(predictor)) %>%
      mean_qi(.epred)
    pred_global_summary$linetype <- pop_credible
    
    # 2) Species-specific prediction grid (restricted to observed range)
    species_pars <- grep(paste0("r_", group, ".*", predictor), colnames(post), value = TRUE)
    
    species_credible <- sapply(species_pars, function(sp) {
      sp_post <- post[[paste0("b_", predictor)]] + post[[sp]]  # total slope
      if(quantile(sp_post, 0.025) <= 0 & quantile(sp_post, 0.975) >= 0) "dashed" else "solid"
    })
    
    species_linetypes <- data.frame(species = species_names, linetype = species_credible[2:6])
    
    species_ranges <- data %>%
      group_by(.data[[group]]) %>%
      summarise(
        min_pred = min(.data[[predictor]]),
        max_pred = max(.data[[predictor]]),
        .groups = "drop"
      )
    
    pred_species_grid <- species_ranges %>%
      rowwise() %>%
      mutate(
        !!predictor := list(seq(min_pred, max_pred, length.out = 50))
      ) %>%
      unnest(cols = !!predictor)
    
    pred_species <- add_epred_draws(
      pred_species_grid,
      object = fit,
      re_formula = NULL
    )
    
    pred_species_summary <- pred_species %>%
      group_by(.data[[group]], !!sym(predictor)) %>%
      mean_qi(.epred) %>%
      left_join(species_linetypes, by = c("species"))
    
    # Plot with global line + CI and species lines
    p <- ggplot(data, aes_string(x = predictor, y = response)) +
      geom_point(alpha = 0.55, color = 'gray80') +
      
      # # Population-level 95% CI ribbon
      # stat_lineribbon(data = pred_global_summary,
      #             aes_string(x = predictor, ymin = ".lower", ymax = ".upper"),
      #             fill = "grey70", alpha = 0.4) +
      # Population-level line
      geom_line(data = pred_global_summary,
                aes_string(x = predictor, y = ".epred", linetype = 'linetype'),
                linewidth = 1.2, color = "black") +
      # Species lines
      geom_line(data = pred_species_summary,
                aes_string(x = predictor, y = ".epred", linetype = "linetype"),
                color = 'gray80',
                linewidth = 1) +
      
      # scale_color_manual(values = c("#42049EFF","#8204A7FF","#B6308BFF","#F79143FF","#FCCE25FF"), labels = species_labels) +
      scale_linetype_identity(guide = 'none') +
      theme_classic()
    
    p
    
    return(list(
      fit = fit,
      key_summary = key_summary,
      slope_table = slope_table,
      plot = p
    ))
  }
  
  # Loop over all trait pairs
  results <- map(trait_list, fit_single_pair)
  names(results) <- sapply(trait_list, function(x) paste(x[2], "vs", x[1], sep = "_"))
  
  return(results)
}






