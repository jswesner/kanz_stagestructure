source("code/functions.r") # functions
source("code/make_data.r") # make data and load packages 

brm_all <- readRDS(file = "models/species_models/brm_all.rds")
brm_chiro_stages <- readRDS(file = "models/species_models/brm_chiro_stages.rds")


# refit with wider priors

brm_all_wider <- update(brm_all, prior = c(prior(normal(0, 10), class = "Intercept"),
                                           prior(normal(-1, 2), class = "b"),
                                           prior(gamma(1, 0.1), class = "shape"),
                                           prior(exponential(2), class = "sd")),
                        newdata = brm_all$data,
                        chains = 1, iter = 1000,
                        file = "models/brm_all_wider.rds",
                        file_refit = "on_change"
                        )

brm_chiro_wider <- update(brm_chiro_stages, prior = c(prior(normal(0, 10), class = "Intercept"),
                                           prior(normal(0, 4), class = "b"),
                                           prior(gamma(1, 0.1), class = "shape"),
                                           prior(exponential(2), class = "sd")),
                        newdata = brm_chiro_stages$data,
                        chains = 1, iter = 1000,
                        file = "models/brm_chiro_wider.rds",
                        file_refit = "on_change"
)

# compare


# plots of total prey
post_all_a <- brm_all$data %>% 
  distinct(prey_feeding, prey_ecosystem, fish_species) %>% 
  mutate(site = NA, date = NA) %>% 
  add_epred_draws(brm_all, re_formula = NULL, dpar = T) %>% 
  ungroup() %>% 
  select(-site, -date, -.row, -.chain, -.iteration,  -mu, -hu, -shape) %>% 
  mutate(model = "main all prey")

post_all_wider <- brm_all_wider$data %>% 
  distinct(prey_feeding, prey_ecosystem, fish_species) %>% 
  mutate(site = NA, date = NA) %>% 
  add_epred_draws(brm_all_wider, re_formula = NULL, dpar = T) %>% 
  ungroup() %>% 
  select(-site, -date, -.row, -.chain, -.iteration,  -mu, -hu, -shape) %>% 
  mutate(model = "wider all prey")


post_chiro_a <- brm_chiro_stages$data %>% 
  distinct(prey_stage, fish_species) %>% 
  mutate(site = NA, date = NA) %>% 
  add_epred_draws(brm_chiro_stages, re_formula = NULL, dpar = T) %>% 
  ungroup() %>% 
  select(-site, -date, -.row, -.chain, -.iteration,  -mu, -hu, -shape) %>% 
  mutate(model = "main chiro") %>% 
  mutate(prey_feeding = case_when(prey_stage == "l" ~ "consumer",
                                  TRUE ~ "non_consumer"))


post_chiro_wider <- brm_chiro_stages$data %>% 
  distinct(prey_stage, fish_species) %>% 
  mutate(site = NA, date = NA) %>% 
  add_epred_draws(brm_chiro_stages, re_formula = NULL, dpar = T) %>% 
  ungroup() %>% 
  select(-site, -date, -.row, -.chain, -.iteration,  -mu, -hu, -shape) %>% 
  mutate(model = "wider chiro") %>% 
  mutate(prey_feeding = case_when(prey_stage == "l" ~ "consumer",
                                  TRUE ~ "non_consumer"))

# combine

all_combine <- bind_rows(post_all_a, post_all_wider) %>% 
  unite("feed_ecosystem", prey_feeding:prey_ecosystem) %>% 
  group_by(fish_species, .draw, model) %>% 
  mutate(total = sum(.epred)) %>% 
  pivot_wider(names_from = feed_ecosystem, values_from = .epred) %>% 
  mutate(prop_noncons = (non_consumer_aquatic + non_consumer_terrestrial)/total,
         total_aquatic = consumer_aquatic + non_consumer_aquatic,
         prop_noncons_aquatic = non_consumer_aquatic/total_aquatic,
         prop_terrestrial = non_consumer_terrestrial/total) %>% 
  group_by(fish_species) %>% 
  mutate(median = median(prop_noncons)) %>% 
  pivot_longer(cols = c(prop_noncons, prop_noncons_aquatic, prop_terrestrial)) %>% 
  separate(fish_species, c("genus", "species"), remove = F) %>% 
  mutate(genus = str_sub(genus, 1, 1),
         genus_species = paste0(genus, ". ", species),
         genus_species = case_when(fish_species == "Lepisosteidae" ~ "Lepisosteidae",
                                   TRUE ~ genus_species))

prior_sens_plot <- all_combine %>% 
  filter(name == "prop_noncons") %>% 
  mutate(Model = case_when(model == "main all prey" ~ "Main",
                           TRUE ~ "Wider Priors")) %>% 
  ggplot(aes(x = reorder(genus_species, -median), y = value)) + 
  geom_boxplot(aes(fill = Model), outlier.shape = NA) +
  theme_default() +
  labs(y = "Proportion of non-consumer prey",
       x = "Fish taxon") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_grey()

ggsave(prior_sens_plot, file = "plots/prior_sens_plot.jpg", dpi = 500, width = 6.5,
       height = 4)

