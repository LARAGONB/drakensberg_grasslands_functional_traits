################################################################################
# Exploratory Data Analysis of Functional Traits from the Drakensberg Rooties' project
################################################################################
#
# Lina Aragón
# linamaragonb@gmail.com
# September 23, 2025
#
# Description
# This script performs an exploratory data analysis of functional traits data from the Drakensberg Rooties' project.
# It includes loading necessary libraries, creating directories, loading and reshaping data, visualizing data
# through various plots, identifying and removing outliers, recalculating certain traits, and saving the cleaned data.
# Additionally, it generates summary tables for the cleaned data.
################################################################################

# 1. Load libraries ----

# install.packages("devtools", "EVR628tools", "tidyverse", "ggplot2", "plotly", "osfr") #install if needed
# devtools::install_github("jcvdav/EVR628tools")
pkgs <- c("devtools", "EVR628tools", "tidyverse", "ggplot2", "plotly", "osfr")
lapply(pkgs, library, character.only = TRUE)
remove(pkgs)

# 2. Create directories ----
#create_dirs() # (This should only be done once by the owner of the repository)

# 3. Download and load data ----
# retrieve raw data files from the OSF project page
osf_retrieve_node('hk2cy') %>%
  osf_ls_files(path = 'v_root_traits/') %>%
  osf_download(path = 'data/raw/', conflicts = 'overwrite')
  
# load data
data_RFT_raw <- read_csv("data/raw/v_PFCT7_clean_root_traits_2023.csv")

# 4. Create wide table ----
data_RFT_raw_wide <- data_RFT_raw |> 
  select(!unit) |> 
  pivot_wider(names_from = traits, values_from = value) #using traits as names and values as values
data_RFT_wide 
write_csv(data_RFT_raw_wide, "data/raw/v_PFCT7_clean_root_traits_2023_wide.csv")

# 5. Explore data ----

## Graphs o changes in FT along the elevational gradient by species
p <- data_RFT_raw |> 
  filter(traits %in% c("rd","bi","srl","rtd","rdmc","sla","ldmc","bgb_agb")) |> 
  ggplot(aes(x = as.factor(elevation_m_asl), y = value, color = species)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
    alpha = 0.8
  ) +
  facet_wrap(~ traits, scales = "free_y") 

ggplotly(p)

## Graphs of overall differences in FT by species
data_RFT_raw |> 
  filter(traits %in% c("rd","bi","srl","rtd","rdmc","sla","ldmc","bgb_agb")) |> 
  ggplot(aes(x = species, y = value, color = species)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
    alpha = 0.8
  ) +
  facet_wrap(~ traits, scales = "free_y") 

## Graphs o changes in FT along the elevational gradient by species w/t outilers?
data_RFT_raw |> 
  filter(traits %in% c("rd","bi","srl","rtd","rdmc","sla","ldmc","bgb_agb")) |> 
  filter(!id %in% c("FBD8952","FCE1581", "FEK5954", "FFO5284")) |> 
  ggplot(aes(x = as.factor(elevation_m_asl), y = value, color = species)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
    alpha = 0.8
  ) +
  facet_wrap(~ traits, scales = "free_y") 

## Graphs of overall differences in FT by species w/t outilers?
data_RFT_raw |> 
  filter(traits %in% c("rd","bi","srl","rtd","rdmc","sla","ldmc","bgb_agb")) |> 
  filter(!id %in% c("FBD8952","FCE1581", "FEK5954", "FFO5284")) |> 
  ggplot(aes(x = species, y = value, color = species)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
    alpha = 0.8
  ) +
  facet_wrap(~ traits, scales = "free_y") 

## Density plots per elevation
data_RFT_raw %>%
  filter(traits %in% c("rd", "bi", "srl", "rtd", "rdmc", "sla", "ldmc", "bgb_agb")) %>%
  ggplot(aes(x = value, fill = species)) +
  geom_density(alpha = 0.6) +
  facet_wrap(~ elevation_m_asl + traits, scales = "free", ncol = 8, nrow = 4) +
  theme_bw() 

## Density plots per trait
data_RFT_raw %>%
  filter(traits %in% c("rd", "bi", "srl", "rtd", "rdmc", "sla", "ldmc", "bgb_agb")) %>%
  ggplot(aes(x = value, fill = species)) +
  geom_density(alpha = 0.6) +
  facet_wrap(~ traits, scales = "free") +
  theme_bw() 


## Scatter plots to explore relationships between traits and find outliers
p1 <- ggplot(data_RFT_raw_wide, aes(x = leaf_wet_mass, y = leaf_dry_mass)) +
  geom_point() 
ggplotly(p1)

p2 <- ggplot(data_RFT_raw_wide, aes(x = aboveground_biomass, y = belowground_biomass)) +
  geom_point() #Scatter plots to explore relationships between traits and find outliers
ggplotly(p2)

#6. Remove outliers, update ldm, recalculate SLA and LDMC, and create growth form ----
#Outliers identified visually
#FF05284 has a very high bgb_agb value due to a very low agb value < bgb
#FEK5954 has a very high rtd value because we could only get one fine root to be analyzed
#FCE1581 has a very high ldmc value due to a typo in the leaf_dry_mass value
#FBD8952 has a very high ldmc value due to a typo in the leaf_dry_mass value


data_RFT_wide <- data_RFT_raw_wide |> 
  filter(!id %in% c("FFO5284","FEK5954")) |> 
  mutate(
    leaf_dry_mass = case_when(
      id == "FCE1581"  ~ 0.033600, 
      id == "FBD8952"  ~ 0.044870,
      TRUE ~ leaf_dry_mass)) |> 
  mutate(
    sla = leaf_area/leaf_dry_mass,
    ldmc = leaf_dry_mass/leaf_wet_mass) |> 
  mutate(
    growth_form = case_when(
      species %in% c("Eragrostis capensis", "Harpochloa falx", "Themeda triandra") ~ "grass",
      species %in% c("Helichrysum pilosellum", "Senecio glaberrimus") ~ "forb",
      TRUE ~ NA),
    family = case_when(
      species %in% c("Eragrostis capensis", "Harpochloa falx", "Themeda triandra") ~ "Poaceae",
      species %in% c("Helichrysum pilosellum", "Senecio glaberrimus") ~ "Asteraceae",
      TRUE ~ NA)) |> 
  relocate(growth_form, .after = species) |> 
  relocate(family, .after = species)
  


data_RFT <- data_RFT_wide |> 
  pivot_longer(cols=10:34,
               names_to = "traits",
               values_to = "value")

#7. Check graphs after changes in dataset ----
p3 <- ggplot(data_RFT_wide_2, aes(x = leaf_wet_mass, y = leaf_dry_mass)) +
  geom_point()
p3


p4 <- data_RFT_2 |> 
  filter(traits %in% c("rd","bi","srl","rtd","rdmc","sla","ldmc","bgb_agb")) |> 
  ggplot(aes(x = as.factor(elevation_m_asl), y = value, color = species)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
    alpha = 0.8
  ) +
  facet_wrap(~ traits, scales = "free_y") #Graphs o changes in FT along the elevational gradient by species

p4
ggplotly(p4)

#8. Save clean data ----
write_csv(data_RFT, "data/processed/v_PFCT7_clean_functional_traits_2023.csv")
write_csv(data_RFT_wide, "data/processed/v_PFCT7_clean_functional_traits_2023_wide.csv")


#7. Create  summary tables for clean data ----

## Summary table per species across all elevations
summary_trait_RFT_spp <- data_RFT |> 
  filter(traits %in% c("rd","bi","srl","rtd","rdmc","sla","ldmc","bgb_agb")) |> 
  group_by(species,traits) |>
  summarise(
    mean = ifelse(all(is.na(value)), NA, mean(value, na.rm = TRUE)),
    sd = ifelse(all(is.na(value)), NA, sd(value, na.rm = TRUE)),
    min = ifelse(all(is.na(value)), NA, min(value, na.rm = TRUE)),
    max = ifelse(all(is.na(value)), NA, max(value, na.rm = TRUE))
  ) 
summary_trait_RFT_spp
write_csv(summary_trait_RFT_spp, "data/output/v_PFCT7_summary_functional_traits_per_spp_2023.csv")

## Summary table per species per elevation
summary_trait_RFT_spp_ele <- data_RFT |> 
  filter(traits %in% c("rd","bi","srl","rtd","rdmc","sla","ldmc","bgb_agb")) |> 
  group_by(species,elevation_m_asl,traits) |>
  summarise(
    mean = ifelse(all(is.na(value)), NA, mean(value, na.rm = TRUE)),
    sd = ifelse(all(is.na(value)), NA, sd(value, na.rm = TRUE)),
    min = ifelse(all(is.na(value)), NA, min(value, na.rm = TRUE)),
    max = ifelse(all(is.na(value)), NA, max(value, na.rm = TRUE))
  ) 
summary_trait_RFT_spp_ele
write_csv(summary_trait_RFT_spp_ele, "data/output/v_PFCT7_summary_functional_traits_per_spp_ele_2023.csv")
