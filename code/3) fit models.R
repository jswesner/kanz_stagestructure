source("code/load_packages.r") # load packages

#load data
guild_diet_multi_drymass <- readRDS("data/guild_diet_multi_drymass.rds")
chiros_only <- readRDS("data/chiros_only.rds")
# load prior data from 2017 REU experiment (Seidel and Wesner, unpublished)
diet_mgdm_2017 <- read_csv("data/prior_diets.csv") 

# Fit models

# all taxa ----------------------------------------------------------------

# prior predictive 
brm_prior <- brm(bf(sample_mg_dm ~ prey_feeding*prey_ecosystem + (1 | site) + (1 | date) + 
                    (1 + prey_feeding*prey_ecosystem | fish_species), 
                  hu ~ 1 + (1 + prey_feeding*prey_ecosystem|fish_species)),
               data = short_prey_terr,
               family = hurdle_gamma(link = "log"),
               prior = c(prior(normal(0, 5), class = "Intercept"),
                         prior(normal(-1, 1), class = "b"),
                         prior(gamma(1, 0.1), class = "shape"),
                         prior(exponential(4), class = "sd")),
               cores = 4,
               sample_prior = "only",
               file_refit = "on_change",
               file = "models/species_models/brm_prior.rds",
               chains = 1, iter = 1000)

# predict mean from fixed and random effects
epred_rands <- brm_prior$data %>% 
  distinct(prey_feeding, prey_ecosystem, site, date, fish_species) %>% 
  add_epred_draws(brm_prior, re_formula = NULL) %>% 
  ungroup() %>% 
  select(-.row, -.chain, -.iteration) %>% 
  group_by(site, date, .draw, fish_species) %>% 
  summarize(total_mg_dm = sum(.epred)) %>% 
  mutate(prediction = "with random effects")

# predict mean from fixed effects only
epred_only <- brm_prior$data %>% 
  distinct(prey_feeding, prey_ecosystem, site, date, fish_species) %>% 
  add_epred_draws(brm_prior, re_formula = NA) %>% 
  ungroup() %>% 
  select(-.row, -.chain, -.iteration) %>% 
  group_by(site, date, .draw, fish_species) %>% 
  summarize(total_mg_dm = sum(.epred)) %>% 
  mutate(prediction = "without random effects")

# simulate individual fish diets
y_preds <- brm_prior$data %>% 
  distinct(prey_feeding, prey_ecosystem, site, date, fish_species) %>% 
  add_predicted_draws(brm_prior, re_formula = NULL) %>% 
  mutate(prediction = "y_pred") %>% 
  ungroup() %>% 
  select(-.row, -.chain, -.iteration) %>% 
  group_by(site, date, .draw, fish_species) %>% 
  summarize(total_mg_dm = sum(.prediction)) %>% 
  mutate(prediction = "y_pred")

# combine all predictions/simulations
all_prior_preds <- bind_rows(epred_rands, epred_only, y_preds)


# plot and compare to previous data
prior_checks <- all_prior_preds %>% 
  ggplot(aes(x = fish_species, y = total_mg_dm, group = interaction(fish_species, date))) +
  geom_violin() +
  # geom_point() +
  facet_wrap(~prediction) +
  scale_y_log10() +
  labs(y = "Total Prey Mass (mgDM per fish)",
       x = "Fish Species",
       subtitle = "Prior predictive distributions") +
  geom_hline(data = diet_mgdm, 
             aes(yintercept = total_mg_dm),
             alpha = 0.2) +
  coord_flip() +
  theme_default()


ggsave(prior_checks, file = "plots/prior_checks.jpg", dpi = 500,
       width = 8.5, height = 4.5)

# fit models --------------------------------------------------------

# load model
brm_all <- readRDS(file = "models/species_models/brm_all.rds")

# fit model
brm_all <- brm(bf(sample_mg_dm ~ prey_feeding*prey_ecosystem + (1 | site) + (1 | date) + 
                    (1 + prey_feeding*prey_ecosystem | fish_species), 
               hu ~ 1 + (1 + prey_feeding*prey_ecosystem|fish_species)),
             data = prey_terr_noncons,
             family = hurdle_gamma(link = "log"),
             prior = c(prior(normal(0, 5), class = "Intercept"),
                       prior(normal(-1, 1), class = "b"),
                       prior(gamma(1, 0.1), class = "shape"),
                       prior(exponential(4), class = "sd")),
             cores = 4,
             file_refit = "on_change",
             file = "models/species_models/brm_all.rds",
             chains = 1, iter = 1000)

brm_all <- update(brm_all, chains = 4, iter = 2000, cores = 4)

saveRDS(brm_all, file = "models/species_models/brm_all.rds")



# chironomid stages -------------------------------------------------------

brm_chiro_stages <- brm(bf(sample_mg_dm ~ prey_stage + (1 + prey_stage|fish_species) + (1|date) + (1|site), 
                           hu ~ 1 + (1 + prey_stage|fish_species)),
                           data = chiros_only,
                           family = hurdle_gamma(link = "log"),
                           prior = c(prior(normal(0, 5), class = "Intercept"),
                                     prior(normal(0, 2), class = "b"),
                                     prior(exponential(2), class = "sd"),
                                     prior(gamma(1, 0.1), class = "shape")),
                           iter = 2000, chains = 4, cores = 4,
                        file = "models/species_models/brm_chiro_stages.rds",
                        file_refit = "on_change")


# Check models ------------------------------------------------------------

pp_all <- pp_check(brm_all, type = "hist") + 
  labs(subtitle = "All Prey") +
  scale_x_continuous(breaks = c(0, 200, 400))

pp_chiro <- pp_check(brm_chiro_stages, type = "hist") + 
  labs(subtitle = "Chironomids Only",
       x = "Prey Mass (mgDM per fish)")

pp_plots <- pp_all/pp_chiro

ggsave(pp_plots, file = "plots/pp_plots.jpg", width = 7, height = 7)



#BY HAND
# y_rep_all <- brm_all$data %>%
#   add_predicted_draws(brm_all, re_formula = NULL, ndraws = 10) %>% 
#   mutate(data = "y_rep") %>% 
#   bind_rows(brm_all$data %>% select(sample_mg_dm) %>% 
#               rename(.prediction = sample_mg_dm) %>% 
#               mutate(.draw = 0,
#                      data = "y"))
# 
# pp_all <- y_rep_all %>% 
#   ggplot(aes(x = .prediction)) + 
#   geom_histogram(aes(group = .draw, fill = data),
#                  bins = 50) +
#   facet_wrap(~.draw) +
#   scale_fill_grey(start = 0, end = 0.8) +
#   theme_default() +
#   labs(x = "Prey Mass (mgDM per fish)",
#        title = "All Prey") +
#   theme(axis.title.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.text.y = element_blank())  +
#   guides(fill = "none")



# y_rep_chiro <- brm_chiro_stages$data %>%
#   add_predicted_draws(brm_chiro_stages, re_formula = NULL, ndraws = 10) %>% 
#   mutate(data = "y_rep") %>% 
#   bind_rows(brm_chiro_stages$data %>% select(sample_mg_dm) %>% 
#               rename(.prediction = sample_mg_dm) %>% 
#               mutate(.draw = 0,
#                      data = "y"))
# 
# pp_chiro <- y_rep_chiro %>% 
#   ggplot(aes(x = .prediction)) + 
#   geom_histogram(aes(group = .draw, fill = data),
#                  bins = 50) +
#   facet_wrap(~.draw) +
#   scale_fill_grey(start = 0, end = 0.8) +
#   theme_default() +
#   labs(x = "Prey Mass (mgDM per fish)",
#        title = "Chironomids only") +
#   theme(axis.title.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.text.y = element_blank()) +
#   guides(fill = "none")
# 
# 
# pp_plots <- pp_all / pp_chiro




