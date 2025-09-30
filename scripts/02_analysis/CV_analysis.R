################################################################################
# PCA analyses
################################################################################
#
# Lina Aragón
# linamaragonb@gmail.com
# September 29, 2025
#
# Description
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

cv_growth_form <- trait_data |> # cv for growth form
  filter(traits %in% c("root_depth", "veg_height", "rd", "bi", "srl", "rtd", "rdmc", "sla", "ldmc", "bgb_agb")) |> 
  group_by(growth_form, traits) |> 
  summarise(cv = (sd(value, na.rm = TRUE) / mean(value, na.rm = TRUE)) * 100, .groups = "drop")

cv_between_spp <- trait_data |> # cv between species
  filter(traits %in% c("root_depth", "veg_height", "rd", "bi", "srl", "rtd", "rdmc", "sla", "ldmc", "bgb_agb")) |> 
  group_by(species, traits) |> 
  summarise(mean_spp = mean(value, na.rm = TRUE), .groups = "drop") |> #traits species means for each species
  group_by(traits) |> 
  summarise(cv = (sd(mean_spp, na.rm = TRUE) / mean(mean_spp, na.rm = TRUE)) * 100) #cv between species

cv_itv_between <- trait_data |> # cv within species along the elevational gradiente (itv between)
  filter(traits %in% c("root_depth", "veg_height", "rd", "bi", "srl", "rtd", "rdmc", "sla", "ldmc", "bgb_agb")) |> 
  group_by(species, elevation_m_asl, traits) |> 
  summarise(mean_spp_ele = mean(value, na.rm = TRUE), .groups = "drop") |>  # traits elevations means for each species
  group_by(species, traits) |> 
  summarise(cv = (sd(mean_spp_ele, na.rm = TRUE) / mean(mean_spp_ele, na.rm = TRUE)) * 100) #cv itv between

cv_itv_within <- trait_data |> # cv within species within each elevational gradient 
  filter(traits %in% c("root_depth", "veg_height", "rd", "bi", "srl", "rtd", "rdmc", "sla", "ldmc", "bgb_agb")) |> 
  group_by(species, elevation_m_asl, traits) |> 
  summarise(cv = (sd(value, na.rm = TRUE) / mean(value, na.rm = TRUE)) * 100)

#Combine all the cv calculation into one table 

cv_long <- bind_rows(cv_growth_form |> mutate(source = "growth_form"),
                     cv_between_spp |> mutate(source = "between_species"),
                     cv_itv_between |> mutate(source = "itv_between"), 
                     cv_itv_within |>  mutate(source = "itv_within")) |> 
  select(traits, cv, source) |> 
  relocate(source, .after = traits) |>
  group_by(traits) |> 
  mutate(
    percent = 100 * cv / sum (cv),
    source = fct_relevel(source, "itv_within", "itv_between", "between_species", "growth_form")
  )

# 4. Graph ----

ggplot(cv_long, aes(x = traits, y = percent, fill = source)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_bar(stat = "identity") +
  scale_y_continuous(
    breaks = seq(0, 100, by = 20),   # breaks every 20 percent: 0, 20, ..., 100
    expand = c(0,0)                  # optional: removes extra space above/below bars
  ) +
  scale_fill_viridis_d(option = "F", direction = -1, begin = 0.1, end = 0.9) +
  labs(
    x = "",
    y = "% of total variance",
    fill = "Source of Variation"
  ) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12))

# 5. Linear mixed model ----

srl_trait_variation <- lmer(srl ~ (1|growth_form/species) + (1|elevation_m_asl), trait_data_wide)
summary(srl_trait_variation)
as_tibble(VarCorr(srl_trait_variation))

lmm_trait_variation <- trait_data |> 
  filter(traits %in% c("root_depth", "rd", "bi", "srl", "rtd", "rdmc", "veg_height", "sla", "ldmc", "bgb_agb" )) |>
  group_by(traits) |>
  nest() |>
  mutate(
    model = map(data, ~ lmer(value ~ (1|growth_form/species) + (1|elevation_m_asl), data = .x)),
    varcomp = map(model, ~ as_tibble(VarCorr(.x)) |> 
                    select(grp,vcov, sdcor) |>
                    mutate(
                      total_var = sum(vcov),
                      proportion = vcov/total_var
                    ))
  )


prop_var <- lmm_trait_variation |> 
  select(varcomp) |> 
  unnest(varcomp)

  

