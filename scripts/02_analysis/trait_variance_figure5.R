################################################################################
# Coefficient of variation and Percentage of total trait variance explained by
# different biologically hierarchichal levels (growth form, species, itv between
# elevation and itv within elevations)
################################################################################
#
# Lina Aragón
# linamaragonb@gmail.com
# September 29, 2025
#
# Description
# In this script we calculated the CV for the biologically hierarchichal level: 
# growth form, species, itv between elevation and itv within elevations.
# Addtionally, we calculated the percentage of total trait variance explained by
# the same biologically hierarchichal levels. For this we used a linear mixed model
################################################################################

# 1. Load libraries ----

#install.packages("devtools", "tidyverse", "ggplot2", "plotly", "vegan", "ggvegan", "ggrepel", "glue", "viridis", "fixest", "lmtest", "corrplot") #install if needed
# devtools::install_github("gavinsimpson/ggvegan")
pkgs <- c("devtools", "tidyverse", "ggplot2", "plotly", "viridis", "fixest", "lmtest", "lme4", "ggh4x")
lapply(pkgs, library, character.only = TRUE)
remove(pkgs)

# 2. Load data ----
trait_data <- read_csv("data/processed/v_PFCT7_clean_functional_traits_2023.csv")
trait_data_wide <- read_csv("data/processed/v_PFCT7_clean_functional_traits_2023_wide.csv")

# 3. Calculate CV = coefficient of variation ---- 
# For growth form, between species, withing species along the elevational gradient, and withing species at each elevation

##cv for all plants ----

cv_all_plants <- trait_data |> 
  filter(traits %in% c("root_depth", "veg_height", "rd", "bi", "srl", "rtd", "rdmc", "sla", "ldmc", "leaf_thickness", "bgb_agb")) |> 
  group_by(traits) |> 
  summarise(cv = (sd(value, na.rm = TRUE) / mean(value, na.rm = TRUE)) * 100, .groups = "drop")

## cv for growth form ----
cv_growth_form <- trait_data |> 
  filter(traits %in% c("root_depth", "veg_height", "rd", "bi", "srl", "rtd", "rdmc", "sla", "ldmc", "leaf_thickness", "bgb_agb")) |> 
  group_by(growth_form, traits) |> 
  summarise(cv = (sd(value, na.rm = TRUE) / mean(value, na.rm = TRUE)) * 100, .groups = "drop")

## cv between species ----
cv_between_spp <- trait_data |> 
  filter(traits %in% c("root_depth", "veg_height", "rd", "bi", "srl", "rtd", "rdmc", "sla", "ldmc", "leaf_thickness", "bgb_agb")) |> 
  group_by(species, traits) |> 
  summarise(mean_spp = mean(value, na.rm = TRUE), .groups = "drop") |> #traits species means for each species
  group_by(traits) |> 
  summarise(cv = (sd(mean_spp, na.rm = TRUE) / mean(mean_spp, na.rm = TRUE)) * 100) #cv between species

## cv within species overall ----
cv_within_spp <- trait_data |>
  filter(traits %in% c("root_depth", "veg_height", "rd", "bi", "srl", "rtd", "rdmc", "sla", "ldmc", "leaf_thickness", "bgb_agb")) |> 
  group_by(species, traits) |>
  summarise(cv = (sd(value, na.rm = TRUE) / mean(value, na.rm = TRUE)) * 100, .groups = "drop")

## cv within species along the elevational gradiente (itv between) ----
cv_itv_between <- trait_data |> 
  filter(traits %in% c("root_depth", "veg_height", "rd", "bi", "srl", "rtd", "rdmc", "sla", "ldmc", "leaf_thickness", "bgb_agb")) |> 
  group_by(species, elevation_m_asl, traits) |> 
  summarise(mean_spp_ele = mean(value, na.rm = TRUE), .groups = "drop") |>  # traits elevations means for each species
  group_by(species, traits) |> 
  summarise(cv = (sd(mean_spp_ele, na.rm = TRUE) / mean(mean_spp_ele, na.rm = TRUE)) * 100) #cv itv between

## cv within species within each elevational gradient ----
cv_itv_within <- trait_data |> 
  filter(traits %in% c("root_depth", "veg_height", "rd", "bi", "srl", "rtd", "rdmc", "sla", "ldmc", "leaf_thickness", "bgb_agb")) |> 
  group_by(species, elevation_m_asl, traits) |> 
  summarise(cv = (sd(value, na.rm = TRUE) / mean(value, na.rm = TRUE)) * 100, .groups = "drop")
  

## Combine all the cv calculation into one table ----

cv_long <- bind_rows(cv_all_plants  |> mutate(source = "all_plants"),
                     cv_growth_form |> mutate(source = "growth_form"),
                     cv_between_spp |> mutate(source = "between_species"),
                     cv_within_spp  |> mutate(source = "within_species"),
                     cv_itv_between |> mutate(source = "itv_between"), 
                     cv_itv_within  |> mutate(source = "itv_within")) |> 
  # select(traits, cv, source) |> 
  relocate(source, .after = traits) |> 
  relocate(growth_form, .after = cv) |> 
  arrange(traits)

### Export table ----
write_csv(cv_long, "results/tab/coefficient_variation_traits.csv")

## Plots for CV ----
### all_within_cv ----
all_within_cv <- cv_long |> 
  filter(source %in% c("all_plants","within_species")) |>
  mutate(
    traits = factor(traits, levels = traits_levels),
    growth_form = case_when(
      species %in% c("Eragrostis capensis", "Harpochloa falx", "Themeda triandra") ~ "Grass",
      species %in% c("Helichrysum pilosellum", "Senecio glaberrimus") ~ "Forb",
      TRUE ~ NA),
    species = if_else(is.na(species), "ALL", species), 
    species = fct_relevel(species, "Themeda triandra", after = 2),
    species = fct_relevel(species, "ALL", after = 5)) |> 
  arrange(traits) |> 
  ggplot(aes(x = species, y = cv, fill = growth_form)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(
    values = c("Forb" = "#6D1F56FF",
               "Grass" = "#B91657FF",
               "NA" = "grey"),
    labels = c("Forb" = "Forb",
               "Grass" = "Grass",
               "NA" = "")) +
  scale_x_discrete(
    labels = c(species_labels, "ALL" = "ALL")) +
  labs(
    y = "Coefficient of Variation (%)",
    x = "") + 
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks = element_line(linewidth = 1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_line(linewidth = 1, colour = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    strip.background = element_rect(fill = "#E0E0E0", color = "black"), # gray background, black border
    strip.text = element_text(color = "black", size = 14)) +
  facet_wrap2(vars(traits), ncol = 3, scales = "free_y", axes = "all", remove_labels = "all",
              labeller = as_labeller(traits_labels))
all_within_cv 

#### Save plot ----
ggsave("results/img/trait_variance/all_within_cv_tiff.tiff", all_within_cv,
       width = 35, height = 20, units = "cm", dpi = 300)
ggsave("results/img/trait_variance/all_within_cv_png.png", all_within_cv,
       width = 35, height = 20, units = "cm", dpi = 300)

### whitin_itvb_itvw ----

# New ordered levels with a placeholder in the 9th position
new_levels <- c(
  "leaf_thickness", "ldmc", "sla",      
  "bi", "rd", "rdmc",                   
  "rtd", "srl",                         
  "EMPTY_PANEL",                        
  "veg_height", "root_depth", "bgb_agb")

# Extend traits_labels so the placeholder has an empty strip label
traits_labels_extended <- traits_labels
traits_labels_extended["EMPTY_PANEL"] <- ""  # no text for the placeholder strip

# Build a fill vector matching the 12 panels.
fills <- c(
  rep(group_colors[1], 3),   
  rep(group_colors[2], 3),   
  rep(group_colors[2], 2),   
  "white",                   
  rep(group_colors[3], 3))
fills_named <- setNames(fills, new_levels)

# Make strip border color mapping and set EMPTY_PANEL strip border to transparent
strip_border_colors <- rep("black", length(new_levels))
names(strip_border_colors) <- new_levels
strip_border_colors["EMPTY_PANEL"] <- "transparent"   # no strip border for empty panel

# build helper df with one row PER FACET (traits column must match facet variable)
border_traits_df <- data.frame(
  traits = setdiff(new_levels, "EMPTY_PANEL"),
  xmin = -Inf, xmax = Inf,
  ymin = -Inf, ymax = Inf,
  stringsAsFactors = FALSE
)

# ensure traits column has same factor levels as in the main data (helps matching)
border_traits_df$traits <- factor(border_traits_df$traits, levels = new_levels)

whitin_itvb_itvw <- cv_long |> 
  filter(source %in% c("within_species", "itv_between", "itv_within")) |>
  mutate(
    traits = factor(traits, levels = new_levels),
    species = fct_relevel(species, "Themeda triandra", after = 2),
    source = factor(source, levels = c("within_species", "itv_between", "itv_within"))) |>
  arrange(traits, source) |>
  ggplot(aes(x = species, y = cv, fill = source)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap2(
    vars(traits),
    ncol = 3,
    scales = "fixed",
    axes = "all",
    remove_labels = "all",
    drop = FALSE,
    labeller = as_labeller(traits_labels_extended),
    strip = strip_themed(
      background_x = elem_list_rect(
        fill = scales::alpha(fills_named, 1),
        color = NA))) +
  scale_fill_manual(
    values = c(
      "within_species" = "gray0",
      "itv_between" = "gray40",
      "itv_within" = "gray80"),
    labels = c(
      "within_species" = "Within species",
      "itv_between" = expression("ITV"["between"]),
      "itv_within" = expression("ITV"["within"]))) +
  scale_x_discrete(
    labels = c(species_labels, "ALL" = "ALL")) +
  labs(
    y = "Coefficient of Variation (%)",
    x = "") + 
  theme_bw(base_size = 16) +
  theme(
    strip.text = element_text(color = "white", size = 14),
    strip.background = element_rect(colour = NA),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks = element_line(colour = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    legend.title = element_blank(),
    legend.text = element_text(size = 16))


whitin_itvb_itvw
#### Save plot ----
ggsave("results/img/trait_variance/whitin_itvb_itvw_cv_png.png", whitin_itvb_itvw,
       width = 35, height = 30, units = "cm", dpi = 300)


# 4. Linear mixed model ---- Where the variation in the trait is found ----

## Table with proportion of variance per growth_form, species, elevation, residual ----
lmm_trait_variation <- trait_data |> 
  filter(traits %in% traits_levels) |>
  group_by(traits) |>
  nest() |>
  mutate(
    model = map(data, ~ lmer(value ~ (1|growth_form) + (1|species) + (1|species:elevation_m_asl), data = .x)),
    varcomp = map(model, ~ as_tibble(VarCorr(.x)) |> 
                    dplyr::select(grp,vcov, sdcor) |>
                    mutate(
                      total_var = sum(vcov),
                      proportion = (vcov/total_var) * 100))) |> 
  dplyr::select(traits, varcomp) |> 
  unnest(varcomp) |> 
  rename(source = grp) |> 
  mutate(
    source = case_when(
      source == "growth_form" ~ "Growth form",
      source == "species" ~ "Species",
      source == "species:elevation_m_asl" ~ "ITV_between",
      source == "Residual" ~ "ITV_within"
    ),
    traits = factor(traits, levels = traits_levels),
    source = factor(source, levels = c("Growth form", "Species", "ITV_between", "ITV_within"))) |>
  arrange(traits, source)

#Comparisons across general (3) hierarchical levels
source_general_model <- lmm_trait_variation |>  
  ungroup() |> 
  mutate(trait_group = factor(traits_groups[as.character(traits)], levels = c("Leaf", "Roots", "Plant size")),
         source_general = case_when(
           source == "ITV_between" ~ "ITV",
           source == "ITV_within" ~ "ITV",
           .default = as.character(source)), .before = source) |> 
  nest(.by = source_general) |>
  mutate(
    model = map(data, \(df) lm(proportion ~ trait_group, data = df)),
    summary = map(model, \(df) summary(df)),
    anova = map(model, \(df) car::Anova(df, type = 3)),
    emmeans = map(model, \(df) emmeans(df, pairwise ~ trait_group, adjust = "bh")))

#Comparisons across all (4) hierarchical levels
source_model <- lmm_trait_variation |>  
  ungroup() |> 
  mutate(trait_group = factor(traits_groups[as.character(traits)], 
                              levels = c("Leaf", "Roots", "Plant size"))) |> 
  nest(.by = source) |>
  mutate(
    model = map(data, \(df) lm(proportion ~ trait_group, data = df)),
    summary = map(model, \(df) summary(df)),
    anova = map(model, \(df) car::Anova(df, type = 3)),
    emmeans = map(model, \(df) emmeans(df, pairwise ~ trait_group, adjust = "bh")))

#Summary across general (3) hierarchical levels
source_summary <- lmm_trait_variation |>  mutate(
  trait_group = factor(traits_groups[as.character(traits)], 
                       levels = c("Leaf", "Roots", "Plant size")),
  source_general = case_when(
    source == "ITV_between" ~ "ITV",
    source == "ITV_within" ~ "ITV",
    .default = as.character(source)), .before = source) |> 
  group_by(trait_group, source_general, source) |> 
  summarise(mean_prop = mean(proportion),
            min_prop = min(proportion),
            max_prop = max(proportion)) |> 
  arrange(trait_group, mean_prop)


#Summary across all (4) hierarchical levels
source_general_summary <- lmm_trait_variation |>  mutate(
  trait_group = factor(traits_groups[as.character(traits)], 
                       levels = c("Leaf", "Roots", "Plant size")),
  source_general = case_when(
    source == "ITV_between" ~ "ITV",
    source == "ITV_within" ~ "ITV",
    .default = as.character(source)), .before = source) |> 
  group_by(trait_group, source_general) |> 
  summarise(mean_prop = mean(proportion),
            min_prop = min(proportion),
            max_prop = max(proportion)) |> 
  arrange(trait_group, mean_prop)

source_detailed_summary <- lmm_trait_variation |>  
  mutate(
    trait_group = factor(traits_groups[as.character(traits)], 
                         levels = c("Leaf", "Roots", "Plant size"))) |>
  group_by(trait_group, source) |>                        # ← source has all 4 levels
  summarise(mean_prop = mean(proportion),
            min_prop  = min(proportion),
            max_prop  = max(proportion),
            .groups = "drop") |> 
  arrange(trait_group, mean_prop)


### Export tables  ----
write_csv(lmm_trait_variation, "results/tab/proportion_total_trait_variance.csv")
write_rds(lmm_trait_variation, "data/output/proportion_total_trait_variance.rds")

write_rds(source_model, "data/output/source_model.rds")  
write_rds(source_general_model, "data/output/source_general_model.rds")

write_csv(source_summary, "results/tab/source_summary_trait_variation.csv")
write_csv(source_general_summary, "results/tab/source_general_summary_trait_variation.csv")
# 5. Visualize model results ----

perc_plot <- ggplot(lmm_trait_variation, aes(x = traits, y = proportion, fill = source)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(
    breaks = seq(0, 100, by = 20),   # breaks every 20 percent: 0, 20, ..., 100
    expand = c(0,0)) +                  # optional: removes extra space above/below bars
  scale_fill_manual(
    values = c(
      "Growth form" = "#F7C5A5FF",
      "Species" = "#ED4F3EFF",
      "ITV_between" = "#931C5BFF",
      "ITV_within" = "#261433FF"),
    labels = c(
      "Growth form" = "Growth form",
      "Species" = "Species",
      "ITV_between" = expression("ITV"["between"]),
      "ITV_within" = expression("ITV"["within"]))) +
  scale_x_discrete(
    labels = traits_labels) +
  labs(
    x = "",
    y = "% of total variance",
    fill = "") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.justification = c("right", "top"))

perc_plot



perc_plot_leaf <- ggplot(lmm_trait_variation |> 
                           filter(traits %in% c("leaf_thickness", "ldmc", "sla")), aes(x = traits, y = proportion, fill = source)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_y_continuous(
    breaks = seq(0, 100, by = 20),   
    expand = c(0,0)) +                
  scale_fill_manual(
    values = group_colors_leaf,
    labels = c(
      "Growth form" = "Growth form",
      "Species" = "Species",
      "ITV_between" = expression("ITV"["between"]),
      "ITV_within" = expression("ITV"["within"]))) +
  scale_x_discrete(
    labels = traits_labels) +
  labs(
    x = "",
    y = "% of total variance",
    fill = "") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.box.spacing = unit(0, "cm"),
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
    legend.title = element_text(angle = 45, vjust = 0.05, hjust = 0.1),
    legend.text = element_text(size = 12),
    legend.justification = c("right", "top"),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

perc_plot_leaf

perc_plot_roots <- ggplot(lmm_trait_variation |> 
                           filter(traits %in% c("bi", "rd", "rdmc", "rtd", "srl")), aes(x = traits, y = proportion, fill = source)) +
  geom_bar(stat = "identity", position = "stack", width = 0.9) +
  scale_y_continuous(
    breaks = seq(0, 100, by = 20),   
    expand = c(0,0)) +                
  scale_fill_manual(
    values = group_colors_roots,
    labels = c(
      "Growth form" = "Growth form",
      "Species" = "Species",
      "ITV_between" = expression("ITV"["between"]),
      "ITV_within" = expression("ITV"["within"]))) +
  scale_x_discrete(
    labels = traits_labels) +
  labs(
    x = "",
    y = "% of total variance",
    fill = "") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.box.spacing = unit(0, "cm"),
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
    legend.title = element_text(angle = 45, vjust = 0.05, hjust = 0.1),
    legend.text = element_text(size = 12),
    legend.justification = c("right", "top"),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

perc_plot_roots

perc_plot_plant_size <- ggplot(lmm_trait_variation |> 
                           filter(traits %in% c("veg_height", "root_depth", "bgb_agb")), aes(x = traits, y = proportion, fill = source)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_y_continuous(
    breaks = seq(0, 100, by = 20),   
    expand = c(0,0)) +                
  scale_fill_manual(
    values = group_colors_plant_size,
    labels = c(
      "Growth form" = "Growth form",
      "Species" = "Species",
      "ITV_between" = expression("ITV"["between"]),
      "ITV_within" = expression("ITV"["within"]))) +
  scale_x_discrete(
    labels = traits_labels) +
  labs(
    x = "",
    y = "% of total variance",
    fill = "") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.box.spacing = unit(0, "cm"),
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
    legend.title = element_text(angle = 45, vjust = 0.05, hjust = 0.1),
    legend.text = element_text(size = 12),
    legend.justification = c("right", "top"),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

perc_plot_plant_size


perc_plot_all_groups <- ggplot(lmm_trait_variation |> 
         mutate(
           trait_group = factor(traits_groups[as.character(traits)], levels = c("Leaf", "Roots", "Plant size")),
           traits = factor(traits, levels = traits_levels),
           source = factor(source, levels = c("Growth form", "Species", "ITV_between", "ITV_within")),
           fill_color = case_when(
             trait_group == "Leaf" & source == "ITV_within" ~ group_colors_leaf["ITV_within"],
             trait_group == "Leaf" & source == "ITV_between" ~ group_colors_leaf["ITV_between"],
             trait_group == "Leaf" & source == "Species" ~ group_colors_leaf["Species"],
             trait_group == "Leaf" & source == "Growth form" ~ group_colors_leaf["Growth form"],
             trait_group == "Roots" & source == "ITV_within" ~ group_colors_roots["ITV_within"],
             trait_group == "Roots" & source == "ITV_between" ~ group_colors_roots["ITV_between"],
             trait_group == "Roots" & source == "Species" ~ group_colors_roots["Species"],
             trait_group == "Roots" & source == "Growth form" ~ group_colors_roots["Growth form"],
             trait_group == "Plant size" & source == "ITV_within" ~ group_colors_plant_size["ITV_within"],
             trait_group == "Plant size" & source == "ITV_between" ~ group_colors_plant_size["ITV_between"],
             trait_group == "Plant size" & source == "Species" ~ group_colors_plant_size["Species"],
             trait_group == "Plant size" & source == "Growth form" ~ group_colors_plant_size["Growth form"],
             TRUE ~ NA_character_)) |> 
         arrange(traits, source), 
       aes(x = traits, y = proportion, fill = fill_color)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_y_continuous(
    limits = c(0, 105),
    breaks = seq(0, 100, by = 20),   
    expand = c(0,0)) +
  scale_fill_identity() +
  scale_x_discrete(
    labels = traits_labels) +
  geom_text(data = cv_all_plants |> mutate(traits = factor(traits, levels = traits_levels)),
            aes(x = traits, y = 102, label = sprintf("%.1f", cv)),
            size = 5, vjust = 0, inherit.aes = FALSE) +
  labs(
    x = "",
    y = "Percentage (%) of total variance",
    fill = "") +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    plot.title = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm"))

perc_plot_all_groups

#Build the legend
color_table <- enframe(c(group_colors_leaf, group_colors_roots, group_colors_plant_size), 
                       name = "source",
                       value = "fill_color") |> 
  mutate(grp = rep(c("Leaf", "Roots", "Plant size"), each = 4),
         grp = factor(grp, levels = c("Leaf", "Roots", "Plant size")),
         source = factor(source, levels = c("ITV_within", "ITV_between", "Species", "Growth form")))

legend_plot <- ggplot(color_table, aes(x = grp, y = source, fill = fill_color)) +
  geom_tile(width = 0.8, height = 0.8) +
  scale_fill_identity() +
  scale_x_discrete(position = "top") +
  scale_y_discrete(position = "right",
                   labels = c(
                     "Growth form" = "Growth form",
                     "Species" = "Species",
                     "ITV_between" = expression("ITV"["between"]),
                     "ITV_within" = expression("ITV"["within"]))) +
  labs(x = "", y = "") +
  coord_fixed(ratio = 1, expand = FALSE) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, size = 14, face = "bold", hjust = 0, vjust = 0,
                               margin = margin(0, 0, 0, 0, unit = "pt")),
    axis.text.y = element_text(size = 14, hjust = 0, margin = margin(0, 0, 0, 0, unit = "pt")),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    plot.title = element_blank(),
    plot.margin = margin(t = 0, r = 0.5, b = 0, l = 0, unit = "pt"))

legend_plot

perc_plot_all_groups_f <- plot_grid(perc_plot_all_groups, legend_plot,
          ncol = 2,
          rel_widths = c(1, 0.3),  # make legend column narrower
          align = "h",
          axis = "tb") +
  theme(panel.border = element_blank(),
        plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(0, 0, 0, 0, unit = "cm"))

perc_plot_all_groups_f

## Save plots ----
ggsave("results/img/trait_variance/perc_plot_all_groups_f_png.png", perc_plot_all_groups_f,
       width = 30, height = 20, units = "cm", dpi = 300)
