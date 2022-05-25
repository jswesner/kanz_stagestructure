library(tidyverse)
library(brms)
library(janitor)
library(readxl)
library(ggridges)
library(lubridate)
library(GGally)
library(patchwork)
library(ggpubr)
library(ggthemes)
library(cowplot)
library(tidybayes)

# load models
brm_all <- readRDS("models/species_models/brm_all.rds")
brm_all_number <- readRDS("models/species_models/brm_all_number.rds")
brm_chiro_stages <- readRDS(file = "models/species_models/brm_chiro_stages.rds")



# Extract Posteriors ------------------------------------------------------

# all
posts_all <- brm_all$data %>% 
  distinct(prey_feeding, prey_ecosystem, fish_species) %>% 
  mutate(site = "na", date = "na") %>% 
  add_epred_draws(brm_all, re_formula = NULL, dpar = T, allow_new_levels = T) %>% 
  ungroup() %>% 
  select(-site, -date, -.row, -.chain, -.iteration,  -mu, -hu, -shape)

# all
posts_number_all <- brm_all_number$data %>% 
  distinct(prey_feeding, prey_ecosystem, fish_species) %>% 
  mutate(site = "na", date = "na") %>% 
  add_epred_draws(brm_all_number, re_formula = NULL, dpar = T, allow_new_levels = T) %>% 
  ungroup() %>% 
  select(-site, -date, -.row, -.chain, -.iteration,  -mu, -hu, -shape)

# chiros only
posts_chiro <- brm_chiro_stages$data %>% 
  distinct(prey_stage, fish_species) %>% 
  mutate(site = "na", date = "na") %>% 
  add_epred_draws(brm_all, re_formula = NULL, dpar = T, allow_new_levels = T) %>% 
  ungroup() %>% 
  select(-site, -date, -.row, -.chain, -.iteration,  -mu, -hu, -shape)

saveRDS(posts_all, file = "posteriors/posts_all.rds")
saveRDS(posts_number_all, file = "posteriors/posts_number_all.rds")
saveRDS(posts_chiro, file = "posteriors/posts_chiro.rds")

# remove large model files
rm(brm_all)
rm(brm_all_number)
rm(brm_chiro_stages)

# post summaries
prop_all_species <- posts_all %>% 
  unite("feed_ecosystem", prey_feeding:prey_ecosystem) %>% 
  group_by(fish_species, .draw) %>% 
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


prop_all_number_species <- posts_number_all %>% 
  unite("feed_ecosystem", prey_feeding:prey_ecosystem) %>% 
  group_by(fish_species, .draw) %>% 
  mutate(total = sum(.epred)) %>% 
  pivot_wider(names_from = feed_ecosystem, values_from = .epred) %>% 
  mutate(prop_noncons = (non_consumer_aquatic + non_consumer_terrestrial)/total,
         total_aquatic = consumer_aquatic + non_consumer_aquatic,
         prop_noncons_aquatic = non_consumer_aquatic/total_aquatic,
         prop_terrestrial = non_consumer_terrestrial/total) %>% 
  group_by(fish_species) %>% 
  pivot_longer(cols = c(prop_noncons, prop_noncons_aquatic, prop_terrestrial)) %>% 
  separate(fish_species, c("genus", "species"), remove = F) %>% 
  mutate(genus = str_sub(genus, 1, 1),
         genus_species = paste0(genus, ". ", species),
         genus_species = case_when(fish_species == "Lepisosteidae" ~ "Lepisosteidae",
                                   TRUE ~ genus_species)) %>% 
  left_join(prop_all_species %>% distinct(fish_species, median))


prop_chiro_species <- posts_chiro %>% 
  group_by(fish_species, .draw) %>% 
  mutate(total = sum(.epred)) %>% 
  pivot_wider(names_from = prey_stage, values_from = .epred) %>% 
  mutate(prop_pup_adult = (p + a)/total) %>% 
  group_by(fish_species) %>% 
  separate(fish_species, c("genus", "species"), remove = F) %>% 
  mutate(genus = str_sub(genus, 1, 1),
         genus_species = paste0(genus, ". ", species),
         genus_species = case_when(fish_species == "Lepisosteidae" ~ "Lepisosteidae",
                                   TRUE ~ genus_species)) %>% 
  mutate(name = "prop_chiro_pup_adult") %>% 
  rename(value = prop_pup_adult) 


prop_all_average <- prop_all_species %>% 
  group_by(.draw, name) %>% 
  select(.draw, name, value) %>% 
  summarize(value = median(value, na.rm = T)) %>% 
  mutate(fish_species = "Average",
         genus_species = "Average",
         median = 0)

prop_all_average_number <- prop_all_number_species %>% 
  group_by(.draw, name) %>% 
  select(.draw, name, value) %>% 
  summarize(value = median(value, na.rm = T)) %>% 
  mutate(fish_species = "Average",
         genus_species = "Average",
         median = 0)

prop_chiro_average <- prop_chiro_species %>% 
  group_by(.draw) %>% 
  select(.draw, value) %>% 
  summarize(value = median(value, na.rm = T)) %>% 
  mutate(fish_species = "Average",
         genus_species = "Average",
         median = 0,
         name = "d) prop_chiro_pup_adult")


prop_posts <- bind_rows(prop_all_species %>% 
                          left_join(prop_all_species %>% distinct(fish_species, median)) , 
                        prop_chiro_species %>% 
                          left_join(prop_all_species %>% distinct(fish_species, median)) ,
                        prop_all_average,
                        prop_chiro_average) %>% 
  mutate(name2 = case_when(grepl("noncons_aquatic", name) ~ "c) Non-consumer prey (aquatic prey only)",
                           grepl("chir", name) ~ "d) Non-consumer prey (chironomids only)",
                           grepl("terres", name) ~ "b) Terrestrial prey (overall)",
                           TRUE ~ "a) Non-consumer prey (overall)")) 


saveRDS(prop_posts, file = "posteriors/prop_posts.rds")

prop_posts_number <- prop_all_number_species  %>%  
  bind_rows(prop_all_average_number) %>% 
  mutate(name2 = case_when(grepl("noncons_aquatic", name) ~ "c) Non-consumer prey (aquatic prey only)",
                           grepl("chir", name) ~ "d) Non-consumer prey (chironomids only)",
                           grepl("terres", name) ~ "b) Terrestrial prey (overall)",
                           TRUE ~ "a) Non-consumer prey (overall)")) 

saveRDS(prop_posts_number, file = "posteriors/prop_posts_number.rds")


# make data to overlay on posteriors

chiros_only <- readRDS("data/chiros_only.rds")

prey_terr_noncons <- guild_diet_multi_drymass %>% 
  group_by(prey_feeding, prey_ecosystem, site, date, fish_species,
           sample_id) %>% 
  summarize(sample_mg_dm = sum(sample_mg_dm, na.rm = T),
            sample_number = sum(number, na.rm = T)) %>% 
  arrange(sample_mg_dm)

prey_all_data <- prey_terr_noncons %>% 
  group_by(fish_species, sample_id, site, date) %>% 
  mutate(total = sum(sample_mg_dm)) %>% 
  unite("feed_ecosystem", c(prey_feeding,prey_ecosystem)) %>% 
  pivot_wider(names_from = feed_ecosystem, values_from = sample_mg_dm) %>% 
  mutate(prop_noncons = (non_consumer_aquatic + non_consumer_terrestrial)/total,
         total_aquatic = consumer_aquatic + non_consumer_aquatic,
         prop_noncons_aquatic = non_consumer_aquatic/total_aquatic,
         prop_terrestrial = non_consumer_terrestrial/total) %>% 
  pivot_longer(cols = c(prop_noncons, prop_noncons_aquatic, prop_terrestrial)) %>% 
  separate(fish_species, c("genus", "species"), remove = F) %>% 
  mutate(genus = str_sub(genus, 1, 1),
         genus_species = paste0(genus, ". ", species),
         genus_species = case_when(fish_species == "Lepisosteidae" ~ "Lepisosteidae",
                                   TRUE ~ genus_species)) %>% 
  left_join(prop_all_species %>% distinct(fish_species, median)) 

prop_data_numeric <- prey_terr_noncons %>% 
  group_by(fish_species, sample_id, site, date) %>% 
  mutate(total = sum(sample_number)) %>% 
  unite("feed_ecosystem", c(prey_feeding,prey_ecosystem)) %>% 
  pivot_wider(names_from = feed_ecosystem, values_from = sample_number) %>% 
  mutate(prop_noncons = (non_consumer_aquatic + non_consumer_terrestrial)/total,
         total_aquatic = consumer_aquatic + non_consumer_aquatic,
         prop_noncons_aquatic = non_consumer_aquatic/total_aquatic,
         prop_terrestrial = non_consumer_terrestrial/total) %>% 
  pivot_longer(cols = c(prop_noncons, prop_noncons_aquatic, prop_terrestrial)) %>% 
  separate(fish_species, c("genus", "species"), remove = F) %>% 
  mutate(genus = str_sub(genus, 1, 1),
         genus_species = paste0(genus, ". ", species),
         genus_species = case_when(fish_species == "Lepisosteidae" ~ "Lepisosteidae",
                                   TRUE ~ genus_species)) %>% 
  left_join(prop_all_species %>% distinct(fish_species, median)) 

chiros_only_data <- chiros_only %>% 
  group_by(fish_species, sample_id, site, date) %>% 
  mutate(total = sum(sample_mg_dm)) %>% 
  select(date, sample_id, site, fish_species, total, prey_stage, sample_mg_dm) %>% 
  pivot_wider(names_from = prey_stage, values_from = sample_mg_dm) %>% 
  mutate(value = (p + a)/total,
         name = "prop_chiro_pup_adult") %>%
  separate(fish_species, c("genus", "species"), remove = F) %>% 
  mutate(genus = str_sub(genus, 1, 1),
         genus_species = paste0(genus, ". ", species),
         genus_species = case_when(fish_species == "Lepisosteidae" ~ "Lepisosteidae",
                                   TRUE ~ genus_species)) %>% 
  left_join(prop_all_species %>% distinct(fish_species, median)) 


prop_data <- bind_rows(prey_all_data, chiros_only_data) %>%
  bind_rows(prey_all_data %>% mutate(genus_species == "Average",
                                     median = 0),
            chiros_only_data %>% mutate(genus_species == "Average",
                                        median = 0)) %>% 
  mutate(name2 = case_when(grepl("noncons_aquatic", name) ~ "c) Non-consumer prey (aquatic prey only)",
                           grepl("chir", name) ~ "d) Non-consumer prey (chironomids only)",
                           grepl("terres", name) ~ "b) Terrestrial prey (overall)",
                           TRUE ~ "a) Non-consumer prey (overall)")) 


saveRDS(prop_data, file = "data/prop_data.rds")
saveRDS(prop_data_numeric, file = "data/prop_data_numeric.rds")
