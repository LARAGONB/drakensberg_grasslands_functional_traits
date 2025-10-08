################################################################################
# Downloading Coordinate Data from the Drakensberg Rooties' project
################################################################################
#
# Lina Aragón - Joe White
# linamaragonb@gmail.com
# October 8, 2025
#
# Description
# This script downloads and processes the coordinate data from the Drakensberg Rooties' project.
################################################################################

# 1. Load libraries ----

# install.packages("tidyverse", "ggplot2", "osfr") #install if needed
pkgs <- c("devtools", "tidyverse", "ggplot2", "osfr", 'sf')
lapply(pkgs, library, character.only = TRUE)
remove(pkgs)

# 2. Download and load data ----
# retrieve coordinate data files from the OSF project page
osf_retrieve_node('hk2cy') %>%
  osf_ls_files(path = '0_coordinates/') %>%
  filter(str_ends(name, '.csv')) %>% 
  osf_download(path = 'data/raw/', conflicts = 'overwrite')

# load data in
coords <- read_csv('data/raw/0_PFCT7_clean_coordinates_2023.csv')
names(coords)

# 3. Summarise data ----
coords_filt <- coords %>%
  filter(aspect == 'west') %>% # only take values from west aspects
  filter(!site_id == 6) %>%
  filter(plot_id == 1)

# 4. Save coordinates as shapefile
write_csv(coords_filt, 'data/processed/site_coordinates.csv')
