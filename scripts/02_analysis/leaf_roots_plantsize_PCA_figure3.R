################################################################################
# PCA analyses
################################################################################
#
# LAB
# September 24, 2025
#
# Description: This code analyses leaf, roots and plant size functional spaces
################################################################################

# SET UP #######################################################################

# Load packages ----------------------------------------------------------------
# install.packages("devtools", "tidyverse", "ggplot2", "plotly", "vegan", "ggvegan", "ggrepel", "glue", "viridis", "fixest", "lmtest", "corrplot") #install if needed
# devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
# remotes::install_github("gavinsimpson/ggvegan")
pkgs <- c("devtools", "tidyverse", "ggplot2", "plotly", "vegan", "ggvegan", "ggrepel", "glue", "viridis",
          "fixest", "lmtest", "corrplot", "FactoMineR", "factoextra", "BiodiversityR",
          "cowplot", "pairwiseAdonis", "car", "tidytext", "rstatix", "emmeans")
# lapply(pkgs, install.packages, character.only = TRUE)
lapply(pkgs, library, character.only = TRUE)
remove(pkgs)


# Load data --------------------------------------------------------------------
trait_data <- read_csv("data/processed/v_PFCT7_clean_functional_traits_2023.csv")
trait_data_wide <- read_csv("data/processed/v_PFCT7_clean_functional_traits_2023_wide.csv")

# PROCESSING ###################################################################
## Roots -----------------------------------------------------------------------
### Selec traits
roots_traits <- trait_data_wide |> 
  select(rd, bi, srl, rtd, rdmc) |> 
  rename_with(toupper)

### Check for skewness and transform traits if needed
roots_sum <- roots_traits |> 
  pivot_longer(everything(), names_to = "trait", values_to = "value") |> 
  group_by(trait) |> 
  summarise(min = min(value),
            max = max(value),
            mean = mean(value),
            median = median(value),
            zero = sum(value == 0),
            skewness = e1071::skewness(value)) ; roots_sum

roots_trans <- roots_traits |> 
  mutate(
    SRL_log = log(SRL),
    RTD_log = log(RTD)) |> 
  select(-SRL, -RTD) |> 
  rename(SRL = SRL_log, RTD = RTD_log)

roots_sum_trans <- roots_trans |> 
  pivot_longer(everything(), names_to = "trait", values_to = "value") |> 
  group_by(trait) |> 
  summarise(min = min(value),
            max = max(value),
            mean = mean(value),
            median = median(value),
            zero = sum(value == 0),
            skewness = e1071::skewness(value)); roots_sum_trans

## Leaf -----------------------------------------------------------------------
### Selec traits
leaf_traits <- trait_data_wide |> 
  select(sla, ldmc, leaf_thickness) |> 
  rename_with(toupper) |> 
  rename(LT = LEAF_THICKNESS)

### Check for skewness and transform traits if needed
leaf_sum <- leaf_traits |> 
  pivot_longer(everything(), names_to = "trait", values_to = "value") |> 
  group_by(trait) |> 
  summarise(min = min(value),
            max = max(value),
            mean = mean(value),
            median = median(value),
            zero = sum(value == 0),
            skewness = e1071::skewness(value)) ; leaf_sum

leaf_trans <- leaf_traits |> 
  mutate(
    LT_log = log(LT),
    SLA_log = log(SLA)) |> 
  select(-LT, -SLA) |> 
  rename(LT = LT_log, SLA = SLA_log)

leaf_sum_trans <- leaf_trans |> 
  pivot_longer(everything(), names_to = "trait", values_to = "value") |> 
  group_by(trait) |> 
  summarise(min = min(value),
            max = max(value),
            mean = mean(value),
            median = median(value),
            zero = sum(value == 0),
            skewness = e1071::skewness(value)); leaf_sum_trans

## Plant size -----------------------------------------------------------------------
### Selec traits
plant_size_traits <- trait_data_wide |> 
  select(root_depth, veg_height, bgb_agb) |> 
  rename(RDepth = root_depth, VHeight = veg_height, BG_AG = bgb_agb)

### Check for skewness and transform traits if needed
plant_size_sum <- plant_size_traits |> 
  pivot_longer(everything(), names_to = "trait", values_to = "value") |> 
  group_by(trait) |> 
  summarise(min = min(value),
            max = max(value),
            mean = mean(value),
            median = median(value),
            zero = sum(value == 0),
            skewness = e1071::skewness(value)) ; plant_size_sum

plant_size_trans <- plant_size_traits |> 
  mutate(
    BG_AG_log = log(BG_AG)) |> 
  select(-BG_AG) |> 
  rename('BG:AG' = BG_AG_log)

plant_size_sum_trans <- plant_size_trans |> 
  pivot_longer(everything(), names_to = "trait", values_to = "value") |> 
  group_by(trait) |> 
  summarise(min = min(value),
            max = max(value),
            mean = mean(value),
            median = median(value),
            zero = sum(value == 0),
            skewness = e1071::skewness(value)); plant_size_sum_trans


# ANALYZE ######################################################################
## Roots -----------------------------------------------------------------------
### PCA ------------------------------------------------------------------------
pca_roots <- roots_trans |> 
  rda(scale = TRUE, center = TRUE)

summary(pca_roots)

### Extract eigenvalues and their percentage
ev_roots <- eigenvals(pca_roots); ev_roots
e_B_roots <- eigenvals(pca_roots)/sum(eigenvals(pca_roots)); e_B_roots
bstick(pca_roots)
screeplot(pca_roots, bstick = T, type = "lines")
### Extract PCs that explain 90% of the variation in the data
k_roots <- which(cumsum(ev_roots) / sum(ev_roots) >= 0.9)[1]; k_roots

### PCA significance
PCAsignificance(pca_roots)
plot1_roots <- ordiplot(pca_roots, choices = c(1,2), scaling = 1)
ordiequilibriumcircle(pca_roots,plot1_roots) #Which traits are more important in each PC

### Tables for ind. loadings and traits loadings
#### Create a table containing the loading of each invidual (sites) in each PC
pca_ind_roots <- as_tibble(bind_cols(
  trait_data_wide |> 
    mutate(species = fct_relevel(species, "Themeda triandra", after = 2)) |> 
    select(id, species, family, growth_form, elevation_m_asl),
  scores(pca_roots, display = "sites", choices = 1:k_roots, scaling = 2)))

#### Create a table containing the loading of each trait (species) in each PC
pca_traits_roots <- scores(pca_roots, display = "species", choices = 1:k_roots, scaling = 2) |> 
  as.data.frame() |> 
  rownames_to_column(var = "traits") 

### PERMANOVA  ----
#### Matrix of distances
roots_distance <- dist(scale(roots_trans), method = "euclidean")

#### Multivarite homogeneity
# elevation
roots_ele_bd <- betadisper(roots_distance, pca_ind_roots$elevation_m_asl)
perm_roots_ele_bd  <- permutest(roots_ele_bd, permutations = 10000)
print(perm_roots_ele_bd)
plot(roots_ele_bd)

# species
roots_spp_bd <- betadisper(roots_distance, pca_ind_roots$species)
perm_roots_spp_bd  <- permutest(roots_spp_bd, permutations = 10000)
print(perm_roots_spp_bd)
plot(roots_spp_bd)

#### Permanova Marginal differences species + elevation
set.seed(4321)
permanova_roots <- adonis2(roots_distance ~ species + elevation_m_asl, 
                           data = as.data.frame(pca_ind_roots), 
                           permutations = 10000,
                           by = "margin")
print(permanova_roots)

spp_per_roots_r2 <- permanova_roots$R2[1] * 100
ele_per_roots_r2 <- permanova_roots$R2[2] * 100
spp_per_roots_p <- format.pval(permanova_roots$`Pr(>F)`[1], eps = 0.05, digits = 3) 
ele_per__p <- format.pval(permanova_roots$`Pr(>F)`[2], eps = 0.05, digits = 3) 


set.seed(4321)
pairwise.adonis2(roots_distance ~ species, 
                 data = as.data.frame(pca_ind_roots),
                 strata = "elevation_m_asl",
                 nperm = 10000,
                 p.adjust.m = "BH")

set.seed(4321)
pairwise.adonis2(roots_distance ~ elevation_m_asl, 
                 data = as.data.frame(pca_ind_roots),
                 strata = "species",
                 nperm = 10000,
                 p.adjust.m = "BH")

### ANOVA PCs ------------------------------------------------------------------
#### squared loadings and percent contribution per axis
pca_traits_roots_load <- pca_traits_roots |> 
  mutate(
    PC1 = PC1^2,
    PC2 = PC2^2,
    PC3 = PC2 ^2,
    CPC1 = PC1/sum(PC1) * 100,
    CPC2 = PC2/sum(PC2) * 100,
    CPC3 = PC3/sum(PC3) * 100) |> 
  pivot_longer(cols = c(CPC1, CPC2, CPC3), names_to = "PC", names_prefix = "C", values_to = "contrib") 

#### Differences among species based on PC scores
model_spp_roots <- pca_ind_roots |> 
  pivot_longer(cols = starts_with("PC"), names_to = "PC", values_to = "scores") |> 
  group_by(PC) |> 
  nest() |> 
  mutate(
    levene = map(data, \(x) levene_test(x, scores ~ species)),
    p_levene = map_dbl(levene, \(x) x$p[1]),
    lm_mod    = map(data, \(x) lm(scores ~ species, data = x)),
    oneway_mod = map(data, \(x) oneway.test(scores ~ species, data = x, var.equal = FALSE)),
    lm_pvalue = map(lm_mod, \(x) car::Anova(x)$`Pr(>F)`[1]),
    oneway_pvalue = map(oneway_mod, \(x) x$p.value),
    W = map(lm_mod, \(x) shapiro.test(residuals(x))$statistic),
    W_pvalue = map(lm_mod, \(x) shapiro.test(residuals(x))$p.value),
    bp_p = map_dbl(lm_mod, \(x) lmtest::bptest(x)$p.value),
    bp_stat = map_dbl(lm_mod, \(x) lmtest::bptest(x)$statistic),
    games_howell = map(data, \(x) games_howell_test(x, scores ~ species)),
    perm_mod = map(data, \(x) lmPerm::lmp(scores ~ species, data = x, perm = "prob")),
    p_perm = map(perm_mod, \(x) summary(x)),
    tukey = map(lm_mod, \(x) emmeans(x, pairwise ~ species, adjust = "BH"))) |> 
  relocate(data, levene, lm_mod, oneway_mod, .after = tukey)

summary(model_spp_roots$lm_mod[[2]])
model_spp_roots$tukey[[2]]

##### Letter for PC1 ~ species plot
letters_pc2_roots <- pca_ind_roots |> 
  group_by(species) |> 
  summarise(y = max(PC2) + 0.25) |> 
  mutate(
    letter = c("ab", "ab", "b", "ab", "a"))

#### Differences among elevation based on PC scores
model_ele_roots <- pca_ind_roots |> 
  pivot_longer(cols = starts_with("PC"), names_to = "PC", values_to = "scores") |> 
  group_by(PC) |> 
  nest() |> 
  mutate(
    lm_mod = map(data, \(x) lm(scores ~ elevation_m_asl, data = x)),
    anova = map(lm_mod, \(x) car::Anova(x)$`Pr(>F)`[1]),
    W = map(lm_mod, \(x) shapiro.test(residuals(x))$statistic),
    W_pvalue = map(lm_mod, \(x) shapiro.test(residuals(x))$p),
    bp_p = map_dbl(lm_mod, \(x) lmtest::bptest(x)$p.value),
    bp_stat = map_dbl(lm_mod, \(x) lmtest::bptest(x)$statistic),
    perm_mod = map(data, \(x) lmPerm::lmp(scores ~ elevation_m_asl, data = x, perm = "prob")),
    p_perm = map_dbl(perm_mod, \(x) coef(summary(x))["elevation_m_asl", "Pr(>|t|)"])) |> 
  relocate(data, lm_mod, perm_mod, .after = p_perm)

## Leaf -----------------------------------------------------------------------
### PCA ------------------------------------------------------------------------
pca_leaf <- leaf_trans |> 
  rda(scale = TRUE, center = TRUE)

summary(pca_leaf)

### Extract eigenvalues and their percentage
ev_leaf <- eigenvals(pca_leaf); ev_leaf
e_B_leaf <- eigenvals(pca_leaf)/sum(eigenvals(pca_leaf)); e_B_leaf
bstick(pca_leaf)
screeplot(pca_leaf, bstick = T, type = "lines")

### Extract PCs that explain 90% of the variation in the data
k_leaf <- which(cumsum(ev_leaf) / sum(ev_leaf) >= 0.9)[1]; k_leaf

### PCA significance
PCAsignificance(pca_leaf)
plot1_leaf <- ordiplot(pca_leaf, choices = c(1,2), scaling = 1)
ordiequilibriumcircle(pca_leaf,plot1_leaf) #Which traits are more important in each PC

### Tables for ind. loadings and traits loadings
#### Create a table containing the loading of each invidual (sites) in each PC
pca_sites_leaf <- as_tibble(bind_cols(
  trait_data_wide |> 
    mutate(species = fct_relevel(species, "Themeda triandra", after = 2)) |> 
    select(id, species, family, growth_form, elevation_m_asl),
  scores(pca_leaf, display = "sites", choices = 1:k_leaf, scaling = 2)))

#### Create a table containing the loading of each trait (species) in each PC
pca_traits_leaf <- scores(pca_leaf, display = "species", choices = 1:k_leaf, scaling = 2) |> 
  as.data.frame() |> 
  rownames_to_column(var = "traits") 

### PERMANOVA  ----
#### Matrix of distances
leaf_distance <- dist(scale(leaf_trans), method = "euclidean")

#### Multivarite homogeneity
# elevation
leaf_ele_bd <- betadisper(leaf_distance, pca_sites_leaf$elevation_m_asl)
perm_leaf_ele_bd  <- permutest(leaf_ele_bd, permutations = 10000)
print(perm_leaf_ele_bd)
plot(leaf_ele_bd)

# species
leaf_spp_bd <- betadisper(leaf_distance, pca_sites_leaf$species)
perm_leaf_spp_bd  <- permutest(leaf_spp_bd, permutations = 10000)
print(perm_leaf_spp_bd)
plot(leaf_spp_bd)

#### Permanova Marginal differences species + elevation
set.seed(4321)
permanova_leaf <- adonis2(leaf_distance ~ species + elevation_m_asl, 
                           data = as.data.frame(pca_sites_leaf), 
                           permutations = 10000,
                           by = "margin")
print(permanova_leaf)

spp_per_leaf_r2 <- permanova_leaf$R2[1] * 100
ele_per_leaf_r2 <- permanova_leaf$R2[2] * 100
spp_per_leaf_p <- format.pval(permanova_leaf$`Pr(>F)`[1], eps = 0.05, digits = 3) 
ele_per_leaf_p <- format.pval(permanova_leaf$`Pr(>F)`[2], eps = 0.05, digits = 3) 


set.seed(4321)
pairwise.adonis2(leaf_distance ~ species, 
                 data = as.data.frame(pca_sites_leaf),
                 strata = "elevation_m_asl",
                 nperm = 10000,
                 p.adjust.m = "BH")

set.seed(4321)
pairwise.adonis2(leaf_distance ~ elevation_m_asl, 
                 data = as.data.frame(pca_sites_leaf),
                 strata = "species",
                 nperm = 10000,
                 p.adjust.m = "BH")

### ANOVA PCs ------------------------------------------------------------------
#### squared loadings and percent contribution per axis
pca_traits_leaf_load <- pca_traits_leaf |> 
  mutate(
    PC1 = PC1^2,
    PC2 = PC2^2,
    CPC1 = PC1/sum(PC1) * 100,
    CPC2 = PC2/sum(PC2) * 100) |> 
  pivot_longer(cols = c(CPC1,CPC2), names_to = "PC", names_prefix = "C", values_to = "contrib") 


#### Differences among species based on PC scores
model_spp_leaf <- pca_sites_leaf |> 
  pivot_longer(cols = starts_with("PC"), names_to = "PC", values_to = "scores") |> 
  group_by(PC) |> 
  nest() |> 
  mutate(
    levene = map(data, \(x) levene_test(x, scores ~ species)),
    p_levene = map_dbl(levene, \(x) x$p[1]),
    lm_mod    = map(data, \(x) lm(scores ~ species, data = x)),
    oneway_mod = map(data, \(x) oneway.test(scores ~ species, data = x, var.equal = FALSE)),
    lm_pvalue = map(lm_mod, \(x) car::Anova(x)$`Pr(>F)`[1]),
    oneway_pvalue = map(oneway_mod, \(x) x$p.value),
    W = map(lm_mod, \(x) shapiro.test(residuals(x))$statistic),
    W_pvalue = map(lm_mod, \(x) shapiro.test(residuals(x))$p.value),
    bp_p = map_dbl(lm_mod, \(x) lmtest::bptest(x)$p.value),
    bp_stat = map_dbl(lm_mod, \(x) lmtest::bptest(x)$statistic),
    games_howell = map(data, \(x) games_howell_test(x, scores ~ species)),
    perm_mod = map(data, \(x) lmPerm::lmp(scores ~ species, data = x, perm = "prob")),
    p_perm = map(perm_mod, \(x) summary(x)),
    tukey = map(lm_mod, \(x) emmeans(x, pairwise ~ species, adjust = "BH"))) |> 
  relocate(data, levene, lm_mod, oneway_mod, .after = tukey)


##### Letter for PC1 ~ species plot
letters_pc1_leaf <- pca_sites_leaf |> 
  group_by(species) |> 
  summarise(y = max(PC1) + 0.08) |> 
  mutate(
    letter = c("a", "b", "c", "d", "e"))

#### Differences among elevation based on PC scores
model_ele_leaf <- pca_sites_leaf |> 
  pivot_longer(cols = starts_with("PC"), names_to = "PC", values_to = "scores") |> 
  group_by(PC) |> 
  nest() |> 
  mutate(
    lm_mod = map(data, \(x) lm(scores ~ elevation_m_asl, data = x)),
    anova = map(lm_mod, \(x) car::Anova(x)$`Pr(>F)`[1]),
    W = map(lm_mod, \(x) shapiro.test(residuals(x))$statistic),
    W_pvalue = map(lm_mod, \(x) shapiro.test(residuals(x))$p),
    bp_p = map_dbl(lm_mod, \(x) lmtest::bptest(x)$p.value),
    bp_stat = map_dbl(lm_mod, \(x) lmtest::bptest(x)$statistic),
    perm_mod = map(data, \(x) lmPerm::lmp(scores ~ elevation_m_asl, data = x, perm = "prob")),
    p_perm = map_dbl(perm_mod, \(x) coef(summary(x))["elevation_m_asl", "Pr(>|t|)"])) |> 
  relocate(data, lm_mod, perm_mod, .after = p_perm)


## Plant_size -----------------------------------------------------------------------
### PCA ------------------------------------------------------------------------
pca_plant_size <- plant_size_trans |> 
  rda(scale = TRUE, center = TRUE)

summary(pca_plant_size)

### Extract eigenvalues and their percentage
ev_plant_size <- eigenvals(pca_plant_size); ev_plant_size
e_B_plant_size <- eigenvals(pca_plant_size)/sum(eigenvals(pca_plant_size)); e_B_plant_size
bstick(pca_plant_size)
screeplot(pca_plant_size, bstick = T, type = "lines")

### Extract PCs that explain 90% of the variation in the data
k_plant_size <- which(cumsum(ev_plant_size) / sum(ev_plant_size) >= 0.9)[1]; k_plant_size

### PCA significance
PCAsignificance(pca_plant_size)
plot1_plant_size <- ordiplot(pca_plant_size, choices = c(1,2), scaling = 1)
ordiequilibriumcircle(pca_plant_size,plot1_plant_size) #Which traits are more important in each PC

### Tables for ind. loadings and traits loadings
#### Create a table containing the loading of each invidual (sites) in each PC
pca_sites_plant_size <- as_tibble(bind_cols(
  trait_data_wide |> 
    mutate(species = fct_relevel(species, "Themeda triandra", after = 2)) |> 
    select(id, species, family, growth_form, elevation_m_asl),
  scores(pca_plant_size, display = "sites", choices = 1:k_plant_size, scaling = 2)))

#### Create a table containing the loading of each trait (species) in each PC
pca_traits_plant_size <- scores(pca_plant_size, display = "species", choices = 1:k_plant_size, scaling = 2) |> 
  as.data.frame() |> 
  rownames_to_column(var = "traits") 

### PERMANOVA  ----
#### Matrix of distances
plant_size_distance <- dist(scale(plant_size_trans), method = "euclidean")

#### Multivarite homogeneity
# elevation
plant_size_ele_bd <- betadisper(plant_size_distance, pca_sites_plant_size$elevation_m_asl)
perm_plant_size_ele_bd  <- permutest(plant_size_ele_bd, permutations = 10000)
print(perm_plant_size_ele_bd)
plot(plant_size_ele_bd)

# species
plant_size_spp_bd <- betadisper(plant_size_distance, pca_sites_plant_size$species)
perm_plant_size_spp_bd  <- permutest(plant_size_spp_bd, permutations = 10000)
print(perm_plant_size_spp_bd)
plot(plant_size_spp_bd)

#### Permanova Marginal differences species + elevation
set.seed(4321)
permanova_plant_size <- adonis2(plant_size_distance ~ species + elevation_m_asl, 
                           data = as.data.frame(pca_sites_plant_size), 
                           permutations = 10000,
                           by = "margin")
print(permanova_plant_size)

set.seed(4321)
pairwise.adonis2(plant_size_distance ~ species, 
                 data = as.data.frame(pca_sites_plant_size),
                 strata = "elevation_m_asl",
                 nperm = 10000,
                 p.adjust.m = "BH")

set.seed(4321)
pairwise.adonis2(plant_size_distance ~ elevation_m_asl, 
                 data = as.data.frame(pca_sites_plant_size),
                 strata = "species",
                 nperm = 10000,
                 p.adjust.m = "BH")

### ANOVA PCs ------------------------------------------------------------------
#### squared loadings and percent contribution per axis
pca_traits_plant_size_load <- pca_traits_plant_size |> 
  mutate(
    PC1 = PC1^2,
    PC2 = PC2^2,
    PC3 = PC2 ^2,
    CPC1 = PC1/sum(PC1) * 100,
    CPC2 = PC2/sum(PC2) * 100,
    CPC3 = PC3/sum(PC3) * 100) |> 
  pivot_longer(cols = c(CPC1, CPC2, CPC3), names_to = "PC", names_prefix = "C", values_to = "contrib") 

#### Differences among species based on PC scores
model_spp_plant_size <- pca_sites_plant_size |> 
  pivot_longer(cols = starts_with("PC"), names_to = "PC", values_to = "scores") |> 
  group_by(PC) |> 
  nest() |> 
  mutate(
    levene = map(data, \(x) levene_test(x, scores ~ species)),
    p_levene = map_dbl(levene, \(x) x$p[1]),
    lm_mod    = map(data, \(x) lm(scores ~ species, data = x)),
    oneway_mod = map(data, \(x) oneway.test(scores ~ species, data = x, var.equal = FALSE)),
    lm_pvalue = map(lm_mod, \(x) car::Anova(x)$`Pr(>F)`[1]),
    oneway_pvalue = map(oneway_mod, \(x) x$p.value),
    W = map(lm_mod, \(x) shapiro.test(residuals(x))$statistic),
    W_pvalue = map(lm_mod, \(x) shapiro.test(residuals(x))$p.value),
    bp_p = map_dbl(lm_mod, \(x) lmtest::bptest(x)$p.value),
    bp_stat = map_dbl(lm_mod, \(x) lmtest::bptest(x)$statistic),
    games_howell = map(data, \(x) games_howell_test(x, scores ~ species)),
    perm_mod = map(data, \(x) lmPerm::lmp(scores ~ species, data = x, perm = "prob")),
    p_perm = map(perm_mod, \(x) summary(x)),
    tukey = map(lm_mod, \(x) emmeans(x, pairwise ~ species, adjust = "BH"))) |> 
  relocate(data, levene, lm_mod, oneway_mod, .after = tukey)

model_spp_plant_size$tukey

##### Letter for PC1 ~ species plot
letters_pc1_plant_size <- pca_sites_plant_size |> 
  group_by(species) |> 
  summarise(y = max(PC1) + 0.1) |> 
  mutate(
    letter = c("b", "ab", "b", "b", "a"))

##### Letter for PC2 ~ species plot
letters_pc2_plant_size <- pca_sites_plant_size |> 
  group_by(species) |> 
  summarise(y = max(PC2) + 0.1) |> 
  mutate(
    letter = c("a", "ab", "c", "d", "bc"))

#### Differences among elevation based on PC scores
model_ele_plant_size <- pca_sites_plant_size |> 
  pivot_longer(cols = starts_with("PC"), names_to = "PC", values_to = "scores") |> 
  group_by(PC) |> 
  nest() |> 
  mutate(
    lm_mod = map(data, \(x) lm(scores ~ elevation_m_asl, data = x)),
    anova = map(lm_mod, \(x) car::Anova(x)$`Pr(>F)`[1]),
    W = map(lm_mod, \(x) shapiro.test(residuals(x))$statistic),
    W_pvalue = map(lm_mod, \(x) shapiro.test(residuals(x))$p),
    bp_p = map_dbl(lm_mod, \(x) lmtest::bptest(x)$p.value),
    bp_stat = map_dbl(lm_mod, \(x) lmtest::bptest(x)$statistic),
    perm_mod = map(data, \(x) lmPerm::lmp(scores ~ elevation_m_asl, data = x, perm = "prob")),
    p_perm = map_dbl(perm_mod, \(x) coef(summary(x))["elevation_m_asl", "Pr(>|t|)"])) |> 
  relocate(data, lm_mod, perm_mod, .after = p_perm)



# VISUALIZE ####################################################################
## Roots -----------------------------------------------------------------------
### PCs contribution ----
contrib_roots <- ggplot(pca_traits_roots_load, aes(x = reorder_within(traits, -contrib, PC), y = contrib, fill = traits)) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 100/ncol(leaf_trans), color = "black", linewidth = 1, linetype = "dashed") +  
  scale_x_reordered() +
  scale_fill_viridis_d(option = "B", begin = 0.1, end = 0.9, alpha = 0.8) +
  facet_wrap(~ PC, scales = "free_x") +  
  labs(x = NULL, 
       y = "% Contribution",
       title = "B. Roots") +
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size = 8),
    axis.ticks = element_line(linewidth = 0.5),
    axis.line = element_line(linewidth = 1, colour = "black"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(),
    legend.key.width = unit(1.2, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.title.position = "plot",
    plot.margin = margin(2, 2, 2, 2, "mm"))

contrib_roots

### PCs 1 and 2 ----
pcs12_roots <- pca_ind_roots |> 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_hline(yintercept = 0, color = "grey90", linewidth = 0.5, linetype = "dashed") +     
  geom_vline(xintercept = 0, color = "grey90", linewidth = 0.5, linetype = "dashed") + 
  geom_point(aes(shape = as.factor(elevation_m_asl)), 
             size = 3, color = "grey70") +
  stat_ellipse(aes(linetype = species), 
               linewidth = 0.8, show.legend = TRUE, color = "grey30") +
  geom_segment(data = pca_traits_roots,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.18, "cm")),
               linewidth = 1.2,
               inherit.aes = FALSE,
               color = group_colors["Roots"]) +
  geom_text_repel(data = pca_traits_roots,
                  aes(x = PC1, y = PC2, label = traits),
                  nudge_x = 0.4,
                  size = 5,
                  show.legend = FALSE,
                  fontface = "bold", 
                  color = group_colors["Roots"]) +
  scale_shape_manual(values = elevation_shapes, name = "Elevation",
                     labels = elevation_labels) +
  scale_linetype_manual(values = species_ellipses, 
                        name = "Species",
                        labels = species_labels) +
  labs(x = glue("PCA1 ({round(e_B_roots[1] * 100, 1)}%)"),
       y = glue("PCA2 ({round(e_B_roots[2] * 100, 1)}%)"),
       title = "B. Roots traits") +
  guides(
    shape = guide_legend(override.aes = list(linetype = "blank", 
                                             colour = "grey70", size = 4)),
    linetype = guide_legend(override.aes = list(shape = NA))) +
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(linewidth = 0.5),
    axis.line = element_line(linewidth = 1, colour = "black"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(),
    legend.key.width = unit(1.2, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.title.position = "plot",
    plot.margin = margin(2, 2, 2, 2, "mm"))

pcs12_roots

### Boxplot PC2 scores ~ species ----
boxplot_roots <- ggplot(aes(y = PC2, x = species, color = species), data = pca_ind_roots) +
  geom_hline(yintercept = 0, color = "grey90", linewidth = 0.5, linetype = "dashed") + 
  geom_point(aes(y = PC2)) +
  geom_boxplot(aes(fill = species), alpha = 0.5) +
  geom_text(data = letters_pc2_roots,
            aes(x = species, y = y, label = letter),
            color = "black", size = 6, fontface = "bold",
            inherit.aes = FALSE) +
  scale_color_manual(values = species_colors, labels = species_labels, name = "Species") +
  scale_fill_manual(values = species_colors, labels = species_labels, , name = "Species") +
  scale_x_discrete(labels = species_labels) +
  labs(
    x = "",
    y = "PC2 scores",
    title = "B. Roots traits") +
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(linewidth = 1),
    axis.line = element_line(linewidth = 1, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.title.position = "plot",
    plot.margin = margin(2,2,2,2, "mm"),
    legend.position = "none")

boxplot_roots

### Regression PC1 scores ~ elevation ----
regression_roots <- ggplot(aes(y = PC1, x = elevation_m_asl), data = pca_ind_roots) +
  geom_point(aes(y = PC1, color = species), size = 3) +
  geom_smooth(method = "lm", color = "black") +
  scale_color_manual(values = species_colors, labels = species_labels, name = "Species") +
  labs(
    x = "",
    y = "PC1 scores",
    title = "B. Roots traits") +
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(linewidth = 1),
    axis.line = element_line(linewidth = 1, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.title.position = "plot",
    plot.margin = margin(2,2,2,2, "mm"))

regression_roots

## Leaf -----------------------------------------------------------------------
### PCs contribution ----
contrib_leaf <- ggplot(pca_traits_leaf_load, aes(x = reorder_within(traits, -contrib, PC), y = contrib, fill = traits)) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 100/ncol(leaf_trans), color = "black", linewidth = 1, linetype = "dashed") +  
  scale_x_reordered() +
  scale_fill_viridis_d(option = "B", begin = 0.1, end = 0.9, alpha = 0.8) +
  facet_wrap(~ PC, scales = "free_x") +  
  labs(x = NULL, 
       y = "% Contribution",
       title = "A. Leaf") +
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size = 8),
    axis.ticks = element_line(linewidth = 0.5),
    axis.line = element_line(linewidth = 1, colour = "black"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(),
    legend.key.width = unit(1.2, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.title.position = "plot",
    plot.margin = margin(2, 2, 2, 2, "mm"))

contrib_leaf

### PCs 1 and 2 ----
pcs12_leaf <- pca_sites_leaf |> 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_hline(yintercept = 0, color = "grey90", linewidth = 0.5, linetype = "dashed") +     
  geom_vline(xintercept = 0, color = "grey90", linewidth = 0.5, linetype = "dashed") + 
  geom_point(aes(shape = as.factor(elevation_m_asl)), 
             size = 3, color = "grey70") +
  stat_ellipse(aes(linetype = species), 
               linewidth = 0.8, show.legend = TRUE, color = "grey30") +
  geom_segment(data = pca_traits_leaf,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.18, "cm")),
               linewidth = 1.2,
               inherit.aes = FALSE,
               color = group_colors["Leaf"]) +
  geom_text_repel(data = pca_traits_leaf,
                  aes(x = PC1, y = PC2, label = traits),
                  nudge_x = 0.4,
                  size = 5,
                  show.legend = FALSE,
                  fontface = "bold", 
                  color = group_colors["Leaf"]) +
  scale_shape_manual(values = elevation_shapes, name = "Elevation",
                     labels = elevation_labels) +
  scale_linetype_manual(values = species_ellipses, 
                        name = "Species",
                        labels = species_labels) +
  labs(x = glue("PCA1 ({round(e_B_leaf[1] * 100, 1)}%)"),
       y = glue("PCA2 ({round(e_B_leaf[2] * 100, 1)}%)"),
       title = "A. Leaf traits") +
  guides(
    shape = guide_legend(override.aes = list(linetype = "blank", 
                                             colour = "grey70", size = 4)),
    linetype = guide_legend(override.aes = list(shape = NA))) +
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(linewidth = 0.5),
    axis.line = element_line(linewidth = 1, colour = "black"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(),
    legend.key.width = unit(1.2, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.title.position = "plot",
    plot.margin = margin(2, 2, 2, 2, "mm"))

pcs12_leaf

### Boxplot PC1 scores ~ species ----
boxplot_leaf <- ggplot(aes(y = PC1, x = species, color = species), data = pca_sites_leaf) +
  geom_hline(yintercept = 0, color = "grey90", linewidth = 0.5, linetype = "dashed") + 
  geom_point(aes(y = PC1)) +
  geom_boxplot(aes(fill = species), alpha = 0.5) +
  geom_text(data = letters_pc1_leaf,
            aes(x = species, y = y, label = letter),
            color = "black", size = 6, fontface = "bold",
            inherit.aes = FALSE) +
  scale_color_manual(values = species_colors, labels = species_labels, name = "Species") +
  scale_fill_manual(values = species_colors, labels = species_labels, , name = "Species") +
  scale_x_discrete(labels = species_labels) +
  labs(
    x = "",
    y = "PC1 scores",
    title = "A. Leaf traits") +
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(linewidth = 1),
    axis.line = element_line(linewidth = 1, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.title.position = "plot",
    plot.margin = margin(2,2,2,2, "mm"),
    legend.position = "none")

boxplot_leaf

### Regression PC2 scores ~ elevation ----
regression_leaf <- ggplot(aes(y = PC2, x = elevation_m_asl), data = pca_sites_leaf) +
  geom_point(aes(y = PC2, color = species), size = 3) +
  geom_smooth(method = "lm", color = "black") +
  scale_color_manual(values = species_colors, labels = species_labels, name = "Species") +
  labs(
    x = "",
    y = "PC2 scores",
    title = "A. Leaf traits") +
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(linewidth = 1),
    axis.line = element_line(linewidth = 1, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.title.position = "plot",
    plot.margin = margin(2,2,2,2, "mm"))

regression_leaf

## Plant_size -----------------------------------------------------------------------
### PCs contribution ----
contrib_plant_size <- ggplot(pca_traits_plant_size_load, aes(x = reorder_within(traits, -contrib, PC), y = contrib, fill = traits)) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 100/ncol(plant_size_trans), color = "black", linewidth = 1, linetype = "dashed") +  
  scale_x_reordered() +
  scale_fill_viridis_d(option = "B", begin = 0.1, end = 0.9, alpha = 0.8) +
  facet_wrap(~ PC, scales = "free_x") +  
  labs(x = NULL, 
       y = "% Contribution",
       title = "C. Plant size") +
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size = 8),
    axis.ticks = element_line(linewidth = 0.5),
    axis.line = element_line(linewidth = 1, colour = "black"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(),
    legend.key.width = unit(1.2, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.title.position = "plot",
    plot.margin = margin(2, 2, 2, 2, "mm"))

contrib_plant_size

### PCs 1 and 2 ----
pcs12_plant_size <- pca_sites_plant_size |> 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_hline(yintercept = 0, color = "grey90", linewidth = 0.5, linetype = "dashed") +     
  geom_vline(xintercept = 0, color = "grey90", linewidth = 0.5, linetype = "dashed") + 
  geom_point(aes(shape = as.factor(elevation_m_asl)), 
             size = 3, color = "grey70") +
  stat_ellipse(aes(linetype = species), 
               linewidth = 0.8, show.legend = TRUE, color = "grey30") +
  geom_segment(data = pca_traits_plant_size,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.18, "cm")),
               linewidth = 1.2,
               inherit.aes = FALSE,
               color = group_colors["Plant size"]) +
  geom_text_repel(data = pca_traits_plant_size,
                  aes(x = PC1, y = PC2, label = traits),
                  nudge_x = 0.4,
                  size = 5,
                  show.legend = FALSE,
                  fontface = "bold", 
                  color = group_colors["Plant size"]) +
  scale_shape_manual(values = elevation_shapes, name = "Elevation",
                     labels = elevation_labels) +
  scale_linetype_manual(values = species_ellipses, 
                        name = "Species",
                        labels = species_labels) +
  labs(x = glue("PCA1 ({round(e_B_plant_size[1] * 100, 1)}%)"),
       y = glue("PCA2 ({round(e_B_plant_size[2] * 100, 1)}%)"),
       title = "C. Plant size traits") +
  guides(
    shape = guide_legend(override.aes = list(linetype = "blank", 
                                             colour = "grey70", size = 4)),
    linetype = guide_legend(override.aes = list(shape = NA))) +
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(linewidth = 0.5),
    axis.line = element_line(linewidth = 1, colour = "black"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(),
    legend.key.width = unit(1.2, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.title.position = "plot",
    plot.margin = margin(2, 2, 2, 2, "mm"))

pcs12_plant_size

### Boxplot PC1 scores ~ species ----
boxplot_pc1_plant_size <- ggplot(aes(y = PC1, x = species, color = species), data = pca_sites_plant_size) +
  geom_hline(yintercept = 0, color = "grey90", linewidth = 0.5, linetype = "dashed") + 
  geom_point(aes(y = PC2)) +
  geom_boxplot(aes(fill = species), alpha = 0.5) +
  geom_text(data = letters_pc1_plant_size,
            aes(x = species, y = y, label = letter),
            color = "black", size = 6, fontface = "bold",
            inherit.aes = FALSE) +
  scale_color_manual(values = species_colors, labels = species_labels, name = "Species") +
  scale_fill_manual(values = species_colors, labels = species_labels, , name = "Species") +
  scale_x_discrete(labels = species_labels) +
  labs(
    x = "",
    y = "PC1 scores",
    title = "C. Plant size traits") +
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(linewidth = 1),
    axis.line = element_line(linewidth = 1, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.title.position = "plot",
    plot.margin = margin(2,2,2,2, "mm"),
    legend.position = "none")

boxplot_pc1_plant_size

### Boxplot PC2 scores ~ species ----
boxplot_pc2_plant_size <- ggplot(aes(y = PC2, x = species, color = species), data = pca_sites_plant_size) +
  geom_hline(yintercept = 0, color = "grey90", linewidth = 0.5, linetype = "dashed") + 
  geom_point(aes(y = PC2)) +
  geom_boxplot(aes(fill = species), alpha = 0.5) +
  geom_text(data = letters_pc2_plant_size,
            aes(x = species, y = y, label = letter),
            color = "black", size = 6, fontface = "bold",
            inherit.aes = FALSE) +
  scale_color_manual(values = species_colors, labels = species_labels, name = "Species") +
  scale_fill_manual(values = species_colors, labels = species_labels, , name = "Species") +
  scale_x_discrete(labels = species_labels) +
  labs(
    x = "",
    y = "PC2 scores",
    title = "D. Plant size traits") +
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(linewidth = 1),
    axis.line = element_line(linewidth = 1, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.title.position = "plot",
    plot.margin = margin(2,2,2,2, "mm"),
    legend.position = "none")

boxplot_pc2_plant_size

### Regression PC1 scores ~ elevation ----
regression_pc1_plant_size <- ggplot(aes(y = PC1, x = elevation_m_asl), data = pca_sites_plant_size) +
  geom_point(aes(y = PC1, color = species), size = 3) +
  geom_smooth(method = "lm", color = "black") +
  scale_color_manual(values = species_colors, labels = species_labels, name = "Species") +
  labs(
    x = "",
    y = "PC1 scores",
    title = "C. Plant size traits") +
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(linewidth = 1),
    axis.line = element_line(linewidth = 1, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.title.position = "plot",
    plot.margin = margin(2,2,2,2, "mm"))

regression_pc1_plant_size

### Regression PC3 scores ~ elevation ----
regression_pc3_plant_size <- ggplot(aes(y = PC3, x = elevation_m_asl), data = pca_sites_plant_size) +
  geom_point(aes(y = PC1, color = species), size = 3) +
  geom_smooth(method = "lm", color = "black") +
  scale_color_manual(values = species_colors, labels = species_labels, name = "Species") +
  labs(
    x = "",
    y = "PC3 scores",
    title = "C. Plant size traits") +
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(linewidth = 1),
    axis.line = element_line(linewidth = 1, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.title.position = "plot",
    plot.margin = margin(2,2,2,2, "mm"))

regression_pc3_plant_size

# EXPORT #######################################################################
## Roots -----------------------------------------------------------------------
ggsave("results/img/pca/pcs12_roots.png", pcs12_roots,
       width = 20, height = 15, units = "cm", dpi = 300)

## Leaf -----------------------------------------------------------------------
ggsave("results/img/pca/pcs12_leaf.png", pcs12_leaf,
       width = 20, height = 15, units = "cm", dpi = 300)

## Plant_size -----------------------------------------------------------------------
ggsave("results/img/pca/pcs12_plant_size.png", pcs12_plant_size,
       width = 20, height = 15, units = "cm", dpi = 300)

## PCA all trait groups
pcs12_roots_nl <- pcs12_roots + theme(legend.position = "none")
pcs12_leaf_nl <- pcs12_leaf + theme(legend.position = "none")
pcs12_plant_size_nl <- pcs12_plant_size + theme(legend.position = "none")

pcs_legend <- get_legend(pcs12_leaf)

pca_all <- plot_grid(pcs12_leaf_nl, 
                     pcs12_roots_nl, 
                     pcs12_plant_size_nl, 
                     pcs_legend) + 
  theme(plot.background = element_rect(fill = "white", color = NA))


ggsave("results/img/pca/pcs12_traits_groups.png", pca_all,
       width = 20, height = 20, units = "cm", dpi = 300)

#Contribution PCs
contrib_all <- plot_grid(contrib_leaf, 
                         contrib_roots, 
                         contrib_plant_size, 
                         ncol = 1) + 
  theme(plot.background = element_rect(fill = "white", color = NA))

ggsave("results/img/pca/contrib_all.png", contrib_all,
       width = 20, height = 20, units = "cm", dpi = 300)

## PCs ~ species all trait groups
pc_spp_all <- plot_grid(boxplot_leaf, boxplot_roots, 
                        boxplot_pc1_plant_size, 
                        boxplot_pc2_plant_size) + 
  theme(plot.background = element_rect(fill = "white", color = NA))

ggsave("results/img/pca/pc_spp_all.png", pc_spp_all,
       width = 20, height = 20, units = "cm", dpi = 300)

## PCs ~ elevation all trait groups
pc_ele_all <- plot_grid(regression_leaf + theme(legend.position = "none"), 
                        regression_roots + theme(legend.position = "none"), 
                        get_legend(regression_leaf),
                        regression_pc1_plant_size + theme(legend.position = "none"), 
                        regression_pc3_plant_size + theme(legend.position = "none"),
                        rel_widths = c(1, 1, 0.25)) + 
  theme(plot.background = element_rect(fill = "white", color = NA))

ggsave("results/img/pca/pc_ele_all.png", pc_ele_all,
      width = 20, height = 20, units = "cm", dpi = 300)


