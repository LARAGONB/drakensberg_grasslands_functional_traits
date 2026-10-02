################################################################################
# Downloading Microclimate Data from the Drakensberg Rooties' project
################################################################################
#
# Lina Aragón - Joe White
# linamaragonb@gmail.com
# October 8, 2025
#
# Description
# This script downloads and processes the site microclimate data from the Drakensberg Rooties' project.
# It includes loading necessary libraries, creating directories, loading and reshaping data, visualizing data
################################################################################

# 1. Load libraries ----

# install.packages("tidyverse", "ggplot2", "osfr") #install if needed
pkgs <- c("devtools", "tidyverse", "ggplot2", "osfr")
lapply(pkgs, library, character.only = TRUE)
remove(pkgs)

# 2. Download and load data ----
# retrieve microclimate data files from the OSF project page
# NOTE: this file is very big 2mil rows, so it is included in gitignore
osf_retrieve_file(id = '68754f967589683d2a5a9b64') %>%
  osf_download(path = 'data/raw/', conflicts = 'overwrite')

# load data in
tomst_microclim <- read_csv('data/raw/xi_PFTC7_clean_elevationgradient_microclimate_2023.csv') %>%
  filter(device == 'Tomst')
names(tomst_microclim)

# 3. Summarise data ----
microclim_sum <- tomst_microclim %>%
  filter(aspect == 'west') %>% # only take values from west aspects
  group_by(elevation_m_asl, climate_variable) %>%
  summarise(mean = mean(value),
            sd = sd(value)) %>%
  mutate(elevation_m_asl = as.factor(elevation_m_asl))

range(tomst_microclim$date_time)[2] - range(tomst_microclim$date_time)[1]

# 4. Visualise data ----
microclim_sum %>%
  ggplot() +
  geom_point(aes(x = elevation_m_asl, y = mean)) +
  facet_grid(~climate_variable)

# 5. Save data ----
write_csv(microclim_sum, 'data/processed/summarised_tomst_microclimate.csv')

