# install.packages('smatr')
library(smatr)
library(tidyverse)

#### Load PCA results for each trait group
pca_all <- read_csv('data/output/PCA_all_traits.csv') %>% rename(allPC1 = PC1, allPC2 = PC2) %>% dplyr::select(id:allPC2)
pca_roots <- read_csv('data/output/PCA_root_traits.csv') %>% rename(rootPC1 = PC1, rootPC2 = PC2) %>% dplyr::select(id:rootPC2)
pca_leaf <- read_csv('data/output/PCA_leaf_traits.csv') %>% rename(leafPC1 = PC1, leafPC2 = PC2) %>% dplyr::select(id:leafPC2)
pca_size <- read_csv('data/output/PCA_size_traits.csv') %>% rename(sizePC1 = PC1, sizePC2 = PC2) %>% dplyr::select(id:sizePC2)

#### organise PCA axes
pca_comb <- pca_all %>%
  left_join(pca_roots, by = names(pca_all)[1:5]) %>%
  left_join(pca_leaf, by = names(pca_all)[1:5]) %>%
  left_join(pca_size, by = names(pca_all)[1:5])

#### run SMA analysis ----
sma_output <- sma(rootPC1 ~ allPC1, data=pca_comb)
# plot root PC1 vs. all PC1
plot(sma_output)

#### create dummy prediction grid
preds <- data.frame(expand.grid(allPC1 = seq(min(pca_comb$allPC1, na.rm = T), max(pca_comb$allPC1, na.rm = T), length.out = 200), stringsAsFactors = FALSE))

# bootstrap data and get predictions
set.seed(1)
preds <- pca_comb %>%
  # create new bootstrapped data sets
  modelr::bootstrap(n = 1000, id = 'boot_num') %>%
  # fit sma to every bootstrap
  group_by(boot_num) %>%
  mutate(., fit = map(strap, ~ sma(rootPC1 ~ allPC1, data=data.frame(.), method="SMA"))) %>%
  ungroup() %>%
  # extract intercept and slope from each fit
  mutate(., intercept = map_dbl(fit, ~coef(.x)[1]),
         slope = map_dbl(fit, ~coef(.x)[2])) %>%
  select(., -fit) %>%
  # get fitted values for each bootstrapped model
  # uses the preds dataframe we made earlier
  group_by(boot_num) %>%
  do(data.frame(fitted = .$intercept + .$slope*preds$allPC1, 
                allPC1 = preds$allPC1)) %>%
  ungroup() %>%
  # calculate the 2.5% and 97.5% quantiles at each allPC1 value
  group_by(., allPC1) %>%
  dplyr::summarise(., conf_low = quantile(fitted, 0.025),
                   conf_high = quantile(fitted, 0.975)) %>%
  ungroup() %>%
  # add fitted value of actual unbootstrapped model
  mutate(., rootPC1 = coef(sma_output)[1] + coef(sma_output)[2]*allPC1)

# extract summary values from the SMA to plot
annotate_summary <- paste0('r = ', ifelse(sma_output$groupsummary$Slope < 0, '-', ''), round(sqrt(sma_output$groupsummary$r2), 2), 
                           if(sma_output$groupsummary$pval < 0.05 & sma_output$groupsummary$pval >= 0.01){
                             '*'
                           } else if(sma_output$groupsummary$pval < 0.01 & sma_output$groupsummary$pval >= 0.001){
                             '**'
                           } else if(sma_output$groupsummary$pval < 0.001){
                             '***'
                           } else{
                             ''
                           })

# plot points, regression, R val
ggplot(pca_comb, aes(allPC1, rootPC1)) +
  geom_point() +
  geom_line(data = preds) +
  geom_ribbon(aes(ymin = conf_low, ymax = conf_high), alpha = 0.1, preds) +
  annotate(geom = 'text', label = annotate_summary, x = max(pca_comb$allPC1)*0.9, y = max(pca_comb$rootPC1)*0.9) +
  theme_bw()

#### Could also run by growth form or species
sma(allPC1 ~ rootPC1*growth_form, data=pca_comb)
plot(sma(allPC1 ~ rootPC1*growth_form, data=pca_comb))

sma(allPC1 ~ rootPC1*species, data=pca_comb)
plot(sma(allPC1 ~ rootPC1*species, data=pca_comb))


