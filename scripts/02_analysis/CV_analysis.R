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
pkgs <- c("devtools", "tidyverse", "ggplot2", "plotly", "viridis", "fixest", "lmtest", "lme4")
lapply(pkgs, library, character.only = TRUE)
remove(pkgs)

# 2. Load data ----
trait_data <- read_csv("data/processed/v_PFCT7_clean_functional_traits_2023.csv")
trait_data_wide <- read_csv("data/processed/v_PFCT7_clean_functional_traits_2023_wide.csv")

# 3. Calculate CV = coefficient of variation ---- 
# For growth form, between species, withing species along the elevational gradient, and withing species at each elevation

##cv for all plants ----

cv_all_plants <- trait_data |> 
  filter(traits %in% c("root_depth", "veg_height", "rd", "bi", "srl", "rtd", "rdmc", "sla", "ldmc", "bgb_agb")) |> 
  group_by(traits) |> 
  summarise(cv = (sd(value, na.rm = TRUE) / mean(value, na.rm = TRUE)) * 100, .groups = "drop")

## cv for growth form ----
cv_growth_form <- trait_data |> 
  filter(traits %in% c("root_depth", "veg_height", "rd", "bi", "srl", "rtd", "rdmc", "sla", "ldmc", "bgb_agb")) |> 
  group_by(growth_form, traits) |> 
  summarise(cv = (sd(value, na.rm = TRUE) / mean(value, na.rm = TRUE)) * 100, .groups = "drop")

## cv between species ----
cv_between_spp <- trait_data |> 
  filter(traits %in% c("root_depth", "veg_height", "rd", "bi", "srl", "rtd", "rdmc", "sla", "ldmc", "bgb_agb")) |> 
  group_by(species, traits) |> 
  summarise(mean_spp = mean(value, na.rm = TRUE), .groups = "drop") |> #traits species means for each species
  group_by(traits) |> 
  summarise(cv = (sd(mean_spp, na.rm = TRUE) / mean(mean_spp, na.rm = TRUE)) * 100) #cv between species

## cv within species overall ----
cv_within_spp <- trait_data |>
  filter(traits %in% c("root_depth", "veg_height", "rd", "bi", "srl", "rtd", "rdmc", "sla", "ldmc", "bgb_agb")) |> 
  group_by(species, traits) |>
  summarise(cv = (sd(value, na.rm = TRUE) / mean(value, na.rm = TRUE)) * 100, .groups = "drop")

## cv within species along the elevational gradiente (itv between) ----
cv_itv_between <- trait_data |> 
  filter(traits %in% c("root_depth", "veg_height", "rd", "bi", "srl", "rtd", "rdmc", "sla", "ldmc", "bgb_agb")) |> 
  group_by(species, elevation_m_asl, traits) |> 
  summarise(mean_spp_ele = mean(value, na.rm = TRUE), .groups = "drop") |>  # traits elevations means for each species
  group_by(species, traits) |> 
  summarise(cv = (sd(mean_spp_ele, na.rm = TRUE) / mean(mean_spp_ele, na.rm = TRUE)) * 100) #cv itv between

## cv within species within each elevational gradient ----
cv_itv_within <- trait_data |> 
  filter(traits %in% c("root_depth", "veg_height", "rd", "bi", "srl", "rtd", "rdmc", "sla", "ldmc", "bgb_agb")) |> 
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
    traits = factor(traits, levels = c("ldmc", "sla", "bi", "rd", "rdmc", "rtd", "srl",
                                       "veg_height","root_depth", "bgb_agb")),
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
    labels = c(
      "Eragrostis capensis" = "ERCA",
      "Harpochloa falx" = "HAFA",
      "Themeda triandra" = "THTR",
      "Helichrysum pilosellum" = "HEPI",
      "Senecio glaberrimus" = "SEGL",
      "ALL" = "ALL")) +
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
  facet_wrap2(vars(traits), ncol = 5, scales = "free_y", axes = "all", remove_labels = "all",
              labeller = as_labeller(c(
                "ldmc" = "LDMC",
                "sla" = "SLA",
                "bi" = "BI",
                "rd" = "RD",
                "rdmc" = "RDMC",
                "rtd" = "RTD",
                "srl" = "SRL",
                "veg_height" = "VHeight",
                "root_depth" = "RDepth",
                "bgb_agb" = "BG:AG")))

#### Save plot ----
ggsave("results/img/all_within_cv_tiff.tiff", all_within_cv,
       width = 40, height = 20, units = "cm", dpi = 300)
ggsave("results/img/all_within_cv_png.png", all_within_cv,
       width = 40, height = 20, units = "cm", dpi = 300)

### whitin_itvb_itvw ----

whitin_itvb_itvw <- cv_long |> 
  filter(source %in% c("within_species", "itv_between", "itv_within")) |>
  mutate(
    traits = factor(traits, levels = c("ldmc", "sla", "bi", "rd", "rdmc", "rtd", "srl",
                                       "veg_height","root_depth", "bgb_agb")),
    species = fct_relevel(species, "Themeda triandra", after = 2),
    source = factor(source, levels = c("within_species", "itv_between", "itv_within"))) |> 
  arrange(traits, source) |> 
  ggplot(aes(x = species, y = cv, fill = source)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(
    values = c(
      "within_species" = "#ED4F3EFF",
      "itv_between" = "#931C5BFF",
      "itv_within" = "#261433FF"),
    labels = c(
      "within_species" = "Within species",
      "itv_between" = expression("ITV"["between"]),
      "itv_within" = expression("ITV"["within"]))) +
  scale_x_discrete(
    labels = c(
      "Eragrostis capensis" = "ERCA",
      "Harpochloa falx" = "HAFA",
      "Themeda triandra" = "THTR",
      "Helichrysum pilosellum" = "HEPI",
      "Senecio glaberrimus" = "SEGL",
      "ALL" = "ALL")) +
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
  facet_wrap2(vars(traits), ncol = 5, scales = "free_y", axes = "all", remove_labels = "all",
              labeller = as_labeller(c(
                "ldmc" = "LDMC",
                "sla" = "SLA",
                "bi" = "BI",
                "rd" = "RD",
                "rdmc" = "RDMC",
                "rtd" = "RTD",
                "srl" = "SRL",
                "veg_height" = "VHeight",
                "root_depth" = "RDepth",
                "bgb_agb" = "BG:AG")))

#### Save plot ----
ggsave("results/img/whitin_itvb_itvw_cv_tiff.tiff", whitin_itvb_itvw,
       width = 40, height = 20, units = "cm", dpi = 300)
ggsave("results/img/whitin_itvb_itvw_cv_png.png", whitin_itvb_itvw,
       width = 40, height = 20, units = "cm", dpi = 300)


# 4. Linear mixed model ---- Where the variation in the trait is found ----

## Table with proportion of variance per growth_form, species, elevation, residual ----
lmm_trait_variation <- trait_data |> 
  filter(traits %in% c("root_depth", "rd", "bi", "srl", "rtd", "rdmc", "veg_height", "sla", "ldmc", "bgb_agb" )) |>
  group_by(traits) |>
  nest() |>
  mutate(
    model = map(data, ~ lmer(value ~ (1|growth_form) + (1|species) + (1|species:elevation_m_asl), data = .x)),
    varcomp = map(model, ~ as_tibble(VarCorr(.x)) |> 
                    select(grp,vcov, sdcor) |>
                    mutate(
                      total_var = sum(vcov),
                      proportion = (vcov/total_var) * 100))) |> 
  select(traits, varcomp) |> 
  unnest(varcomp) |> 
  mutate(
    grp = case_when(
      grp == "growth_form" ~ "Growth form",
      grp == "species" ~ "Species",
      grp == "species:elevation_m_asl" ~ "ITV_between",
      grp == "Residual" ~ "ITV_within"
    ),
    traits = factor(traits, levels = c("ldmc", "sla", "bi", "rd", "rdmc", "rtd", "srl",
                                       "veg_height","root_depth", "bgb_agb")),
    grp = factor(grp, levels = c("Growth form", "Species", "ITV_between", "ITV_within"))) |>
  arrange(traits, grp)

### Export table  ----
write_csv(lmm_trait_variation, "results/tab/proportion_total_trait_variance.csv")

# 5. Visualize model results ----

perc_plot <- ggplot(lmm_trait_variation, aes(x = traits, y = proportion, fill = grp)) +
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
    labels = c(
      "ldmc" = "LDMC",
      "sla" = "SLA",
      "bi" = "BI",
      "rd" = "RD",
      "rdmc" = "RDMC",
      "rtd" = "RTD",
      "srl" = "SRL",
      "veg_height" = "VHeight",
      "root_depth" = "RDepth",
      "bgb_agb" = "BG:AG")) +
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


## Save plots ----

ggsave("results/img/perc_total_variance_tiff.tiff", perc_plot,
       width = 18, height = 10, units = "cm", dpi = 300)
ggsave("results/img/perc_total_variance_png.png", perc_plot,
       width = 18, height = 10, units = "cm", dpi = 300)
