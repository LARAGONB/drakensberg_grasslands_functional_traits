# 1. Load libraries ----

#install.packages("devtools", "tidyverse", "ggplot2", "plotly", "vegan", "ggvegan", "ggrepel", "glue", "viridis", "fixest", "lmtest", "corrplot") #install if needed
# devtools::install_github("gavinsimpson/ggvegan")
pkgs <- c("devtools", "tidyverse", "ggplot2", "plotly", "viridis", "fixest", "lmtest", "lme4", "ggh4x", "emmeans")
lapply(pkgs, library, character.only = TRUE)
remove(pkgs)

# 2. Load data ----
trait_data <- read_csv("data/processed/v_PFCT7_clean_functional_traits_2023.csv")

#Traits levels
traits_levels <- c("leaf_thickness", "ldmc", "sla",
                   "bi", "rd", "rdmc", "rtd", "srl",
                   "veg_height","root_depth", "bgb_agb", 'reproductive_height')

#Traits levels
traits_groups <- c(leaf_thickness = "Leaf", ldmc = "Leaf", sla = "Leaf",
                   bi = "Root", rd = "Root", rdmc = "Root", rtd = "Root", srl = "Root",
                   veg_height = "Plant size", root_depth = "Plant size", bgb_agb = "Plant size", reproductive_height = 'Plant size')

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

# 
# get_var_props <- function(fit) {
#   vc <- as.data.frame(VarCorr(fit))
#   vc$vcov / sum(vc$vcov) * 100
# }
# 
# boot_one_trait <- function(fit, nsim = 1000) {
#   bootMer(fit, FUN = get_var_props, nsim = nsim, 
#           type = "parametric", .progress = "txt")
# }
# 
# # then apply across your nested tibble
# lmm_trait_variation_boot <- trait_data |> 
#   filter(traits %in% traits_levels) |>
#   group_by(traits) |>
#   nest() |>
#   mutate(
#     model = map(data, ~ lmer(value ~ (1|growth_form) + (1|species) + (1|species:elevation_m_asl), data = .x)),
#     boot  = map(model, boot_one_trait),
#     ci    = map(boot, ~ as_tibble(t(apply(.x$t, 2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))))
#   )

###
get_var_props <- function(fit) {
  vc <- as.data.frame(VarCorr(fit))
  props <- vc$vcov / sum(vc$vcov) * 100
  names(props) <- vc$grp
  props
}

boot_one_trait <- function(fit, nsim = 1000) {
  bootMer(fit, FUN = get_var_props, nsim = nsim, 
          type = "parametric", .progress = "txt",
          parallel = "multicore", ncpus = 4)
}

# apply across your nested tibble
lmm_trait_variation_boot <- trait_data |>
  filter(traits %in% traits_levels) |>
  group_by(traits) |>
  nest() |>
  mutate(
    model = map(data, ~ lmer(value ~ (1|growth_form) + (1|species) + (1|species:elevation_m_asl), data = .x)),
    boot  = map(model, boot_one_trait),
    ci    = map(boot, ~ {
      out <- as_tibble(t(apply(.x$t, 2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)),
                       .name_repair = "unique")
      names(out) <- c("lower", "median", "upper")
      out$source <- colnames(.x$t)
      out
    })
  )


# --- Points + whiskers table ---
lmm_ci_points <- lmm_trait_variation_boot |>
  select(traits, ci) |>
  unnest(ci) |>
  mutate(
    source = case_when(
      source == "growth_form" ~ "Growth form",
      source == "species" ~ "Species",
      source == "species:elevation_m_asl" ~ "ITV_between",
      source == "Residual" ~ "ITV_within"
    ),
    traits = factor(traits, levels = traits_levels),
    source = factor(source, levels = c("Growth form", "Species", "ITV_between", "ITV_within"))
  ) |>
  mutate(
    trait_group = traits_groups[as.character(traits)],
    trait_group = factor(trait_group, levels = c("Leaf", "Root", "Plant size"))
  )

# Order traits by group first, then original within-group order
trait_order <- names(traits_groups)[order(match(traits_groups, c("Leaf", "Root", "Plant size")))]

lmm_ci_points <- lmm_ci_points |>
  mutate(traits = factor(traits, levels = trait_order))

####
p_points <- ggplot(lmm_ci_points, aes(x = source, y = median, color = source)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), fatten = 2, linewidth = 0.7, size = 2) +
  facet_nested_wrap(~ trait_group + traits, nrow = 2, nest_line = TRUE) +
  coord_flip() +
  scale_color_manual(values = c(
    "Growth form" = "#1B9E77",
    "Species"     = "#D95F02",
    "ITV_between" = "#7570B3",
    "ITV_within"  = "#E7298A"
  ), guide = "none") +
  labs(y = "Percentage (%) of total variance", x = NULL) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())

p_points

ggsave('results/img/trait_variance/variance_bootstrap.png',
       dpi = 320,
       width = 11.00,
       height = 5.96)
