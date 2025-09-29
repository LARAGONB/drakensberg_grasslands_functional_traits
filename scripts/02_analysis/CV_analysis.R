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

cv_growth_form <- trait_data |> 
  filter(traits %in% c("root_depth", "veg_height", "rd", "bi", "srl", "rtd", "rdmc", "sla", "ldmc", "bgb_agb")) |> 
  group_by(growth_form, traits) |> 
  summarise(cv = (sd(value, na.rm = TRUE) / mean(value, na.rm = TRUE)) * 100, .groups = "drop")

cv_between_spp <- trait_data |> 
  filter(traits %in% c("root_depth", "veg_height", "rd", "bi", "srl", "rtd", "rdmc", "sla", "ldmc", "bgb_agb")) |> 
  group_by(species, traits) |> 
  summarise(mean_spp = mean(value, na.rm = TRUE), .groups = "drop") |> 
  group_by(traits) |> 
  summarise(cv = (sd(mean_spp, na.rm = TRUE) / mean(mean_spp, na.rm = TRUE)) * 100)

cv_itv_between <- trait_data |> 
  filter(traits %in% c("root_depth", "veg_height", "rd", "bi", "srl", "rtd", "rdmc", "sla", "ldmc", "bgb_agb")) |> 
  group_by(species, elevation_m_asl, traits) |> 
  summarise(mean_spp_ele = mean(value, na.rm = TRUE), .groups = "drop") |> 
  group_by(species, traits) |> 
  summarise(cv = (sd(mean_spp_ele, na.rm = TRUE) / mean(mean_spp_ele, na.rm = TRUE)) * 100)

cv_itv_within <- trait_data |> 
  filter(traits %in% c("root_depth", "veg_height", "rd", "bi", "srl", "rtd", "rdmc", "sla", "ldmc", "bgb_agb")) |> 
  group_by(species, elevation_m_asl, traits) |> 
  summarise(cv = (sd(value, na.rm = TRUE) / mean(value, na.rm = TRUE)) * 100)

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


