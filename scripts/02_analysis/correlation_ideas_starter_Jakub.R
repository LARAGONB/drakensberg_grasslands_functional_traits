# Jakub 
# Starter on correlations
# 19/11/2025

# Correlation tests like spearman's p or SMA woud be most appropriate but I don't think they can account for species structure
# Potentially try out Bayesian Gaussian Copula Correlation which itnernet says allows usign species as random effects

# Libraries ----
library(mgcv)
library(ggeffects)


# Data ----
rdata <- trait_data_wide %>% # root traits
  dplyr::select(species, rd, bi, srl, rtd) 

ldata <- trait_data_wide %>% # leaf traits
  dplyr::select(species, sla, ldmc)


# Modelling - not correct, but I was testing out how these analyses work in comparison to your scripts Lina ----
# Gam offer flexibility but also harder to interpret
m_gam2 <- gam(sla ~ s(ldmc) + s(species, bs = "re"), data = ldata)
summary(m_gam2)
plot(m_gam2, pages = 1)
ggeffects::ggpredict(m_gam2, terms = "ldmc [all]") %>% plot()


# We predict SRL and RTD are negatively correlated
# Issue with model not converging with the species as random effect, but species seem to follow independently the general pattern
m_srl_rtd_re <- lmer(srl ~ rtd + (1|species), data = rdata) # not converged
m_srl_rtd <- lm(srl ~ rtd + species, data = rdata)
summary(m_srl_rtd)
ggpredict(m_srl_rtd, terms = "rtd [all]") %>% plot()
ggpredict(m_srl_rtd, terms = "species") %>% plot()

m_gam_srl_rtd <- gam(srl ~ s(rtd) + s(species, bs = "re"), data = rdata)
summary(m_gam_srl_rtd)
plot(m_gam_srl_rtd, pages = 1)
ggeffects::ggpredict(m_gam_srl_rtd, terms = "rtd [all]") %>% plot()

# Correlation test woud be much more appropriate
cor.test(rdata$srl, rdata$rtd, method = "spearman") # -0.75


# SLA and LDMC - example of why not to omit species in the models (relationship reversed!)
m_sla_ldmc_re <- lmer(sla ~ ldmc + (1|species), data = ldata)
summary(m_sla_ldmc_re)
ggeffects::ggpredict(m_sla_ldmc_re, terms = "ldmc [all]") %>% plot()

m_ldmc_sla_re <- lmer(ldmc ~ sla + (1|species), data = ldata)
summary(m_ldmc_sla_re)
ggeffects::ggpredict(m_ldmc_sla_re, terms = "sla [all]") %>% plot()

m_sla_ldmc <- lm(sla ~ ldmc, data = ldata)
summary(m_sla_ldmc)
ggeffects::ggpredict(m_sla_ldmc, terms = "ldmc [all]") %>% plot()

# Here, Spearman's correlation test shows it's not enough
cor.test(ldata$sla, ldata$ldmc, method = "spearman") # 0.61 - positive, while with species it would be negative


# Also, again, if prediction is the goal then GAM offers flexibility but also harder to interpret
m_gam2 <- gam(sla ~ s(ldmc) + s(species, bs = "re"), data = ldata)
summary(m_gam2)
plot(m_gam2, pages = 1)
ggeffects::ggpredict(m_gam2, terms = "ldmc [all]") %>% plot()


