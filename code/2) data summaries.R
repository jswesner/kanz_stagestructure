source("code/load_packages.r") # load packages

guild_diet_multi_drymass <- readRDS("data/guild_diet_multi_drymass.rds")
chiros_only <- readRDS("data/chiros_only.rds")

#fish numbers

fish_table <- guild_diet_multi_drymass %>% distinct(sample_id, fish_species, fish_guild) %>% 
  group_by(fish_species) %>% 
  tally() %>% 
  arrange(-n)

write.csv(fish_table, file = "tables/fish_table.csv", row.names = F)


#proportion of terrestrial prey
prop_terr <- guild_diet_multi_drymass %>% 
  group_by(fish_species, prey_ecosystem, sample_id, site) %>% 
  summarize(total_mg = sum(sample_mg_dm, na.rm = T)) %>% 
  pivot_wider(names_from = prey_ecosystem, values_from = total_mg) %>% 
  mutate(total = sum(aquatic, terrestrial, unknown)) %>% 
  mutate(prop_terrestrial = terrestrial/total) %>% 
  arrange(desc(prop_terrestrial))

prop_nonfeed_aquatic <- guild_diet_multi_drymass %>% filter(prey_ecosystem == "aquatic") %>% 
  group_by(date, fish_guild, prey_feeding, site, fish_species, sample_id) %>% 
  summarize(sample_mg_dm01 = sum(sample_mg_dm) + 0.01) %>% 
  drop_na(prey_feeding) %>% 
  group_by(fish_species, prey_feeding, sample_id, site) %>% 
  summarize(total_mg = sum(sample_mg_dm01, na.rm = T)) %>% 
  pivot_wider(names_from = prey_feeding, values_from = total_mg) %>% 
  mutate(total = sum(consumer, non_consumer)) %>% 
  mutate(prop_aquatic_nonfeeding = non_consumer/total) 


left_join(prop_terr %>% select(-total), prop_nonfeed_aquatic %>% select(-total)) %>% 
  ggplot(aes(fill = site, x = prop_terrestrial, y = prop_aquatic_nonfeeding)) + 
  geom_point() +
  # scale_x_log10() + 
  # scale_y_log10() +
  geom_smooth(method = "lm")

left_join(prop_terr %>% select(-total), prop_nonfeed_aquatic %>% select(-total)) %>% 
  pivot_longer(cols = c(prop_terrestrial, prop_aquatic_nonfeeding)) %>%
  ggplot(aes(x = reorder(fish_species, value), y = value, color = name)) + 
  geom_point() +
  coord_flip()


#mass in fish guts
guild_diet_multi_drymass %>% 
  group_by(sample_id) %>% 
  summarize(total_mg = sum(sample_mg_dm),
            total_number = sum(number)) %>% 
  pivot_longer(-sample_id) %>% 
  group_by(name) %>% 
  summarize(mean = mean(value),
            median = median(value),
            sd = sd(value),
            min = min(value),
            max = max(value))


guild_diet_multi_drymass %>% 
  group_by(sample_id, fish_guild) %>% 
  summarize(total_mg = sum(sample_mg_dm),
            total_number = sum(number)) %>% 
  pivot_longer(cols = c(-sample_id, -fish_guild)) %>% 
  group_by(name, fish_guild) %>% 
  summarize(mean = mean(value),
            median = median(value),
            sd = sd(value),
            min = min(value),
            max = max(value))


#number of fish a prey item was found in
n_fish <- guild_diet_multi_drymass %>% 
  select(sample_id, prey_class, prey_family, prey_taxon, prey_stage, prey_ecosystem, number, sample_mg_dm) %>% 
  group_by(prey_taxon, sample_id) %>% 
  summarize(number = sum(number)) %>% 
  mutate(prey_presence = case_when(number > 0 ~ 1, TRUE ~ 0)) %>% 
  group_by(prey_taxon) %>% 
  summarize(n = sum(prey_presence)) %>% 
  mutate(n_prop = n/max(guild_diet_multi_drymass$sample_id)) %>% 
  arrange(-n_prop)


n_class <- guild_diet_multi_drymass %>% 
  select(sample_id, prey_class, prey_family, prey_taxon, prey_stage, prey_ecosystem, number, sample_mg_dm) %>% 
  group_by(prey_class, sample_id) %>% 
  summarize(number = sum(number)) %>% 
  mutate(prey_presence = case_when(number > 0 ~ 1, TRUE ~ 0)) %>% 
  group_by(prey_class) %>% 
  summarize(n = sum(prey_presence)) %>% 
  mutate(n_prop = n/max(guild_diet_multi_drymass$sample_id))


#proportions
# proportion of prey taxa
totals <- guild_diet_multi_drymass %>% 
  select(prey_class, prey_family, fish_species, fish_guild, prey_taxon, prey_stage, prey_ecosystem, number, sample_mg_dm) %>%
  # group_by(prey_taxon) %>% 
  summarize(number_total = sum(number),
            mg_dm_total = sum(sample_mg_dm)) %>% 
  mutate(prey_taxon = "Grand Total")

diet_summaries <- guild_diet_multi_drymass %>% 
  select(prey_class, prey_family, fish_species, fish_guild, prey_taxon, prey_stage, prey_ecosystem, number, sample_mg_dm) %>%
  group_by(prey_taxon) %>% 
  summarize(number_total = sum(number),
            mg_dm_total = sum(sample_mg_dm)) %>% 
  mutate(number_prop = number_total/sum(number_total),
         mg_prop = mg_dm_total/sum(mg_dm_total)) %>% 
  arrange(-mg_prop) %>% 
  bind_rows(totals)

write.csv(diet_summaries, file = "tables/diet_summaries.csv", row.names = F)


proportions_by_preytaxa <- guild_diet_multi_drymass %>% 
  select(prey_class, prey_family, fish_species, fish_guild, prey_taxon, prey_stage, prey_ecosystem, number, sample_mg_dm) %>%
  group_by(prey_taxon, fish_guild) %>% 
  summarize(number_total = sum(number),
            mg_dm_total = sum(sample_mg_dm)) %>% 
  group_by(fish_guild) %>% 
  mutate(number_prop = number_total/sum(number_total),
         mg_prop = mg_dm_total/sum(mg_dm_total)) %>% 
  arrange(-number_prop) %>% 
  left_join(n_fish)

# proportions without crayfish
proportions_by_preytaxa_nocrayfish <- guild_diet_multi_drymass %>% 
  select(prey_class, prey_family, fish_species, fish_guild, prey_taxon, prey_stage, prey_ecosystem, number, sample_mg_dm) %>%
  group_by(prey_taxon, fish_guild) %>% 
  filter(prey_taxon != "Crayfish") %>% 
  summarize(number_total = sum(number),
            mg_dm_total = sum(sample_mg_dm)) %>% 
  group_by(fish_guild) %>% 
  mutate(number_prop = number_total/sum(number_total),
         mg_prop = mg_dm_total/sum(mg_dm_total)) %>% 
  arrange(-number_prop) %>% 
  left_join(n_fish)


proportions_by_preytaxa %>% 
  mutate(order = mg_prop) %>% 
  pivot_longer(cols = ends_with("prop")) %>% 
  # filter(name != "n_prop") %>% 
  mutate(name = case_when(name == "mg_prop" ~ "% biomass",
                          name == "n_prop" ~ "Occurrence",
                          TRUE ~ "% abundance"),
         name = fct_relevel(name, "% biomass")) %>% 
ggplot(aes(x = reorder(prey_taxon, order), y = value*100, color = name)) +
  geom_point() +
  scale_color_colorblind() +
  geom_line(aes(group = name)) +
  coord_flip() +
  labs(x = "Prey Taxa (ranked by biomass)",
       y = "% in fish diets") +
  # scale_y_log10() +
  theme_bw() +
  theme(legend.title = element_blank()) +
  NULL


#proportions by fish
proportions_by_fishspecies <- guild_diet_multi_drymass %>% 
  select(prey_class, fish_species, prey_family, prey_taxon, prey_stage, prey_ecosystem, number, sample_mg_dm) %>%
  group_by(prey_taxon, fish_species) %>% 
  summarize(number_total = sum(number),
            mg_dm_total = sum(sample_mg_dm)) %>% 
  group_by(fish_species) %>% 
  mutate(number_prop = number_total/sum(number_total),
         mg_prop = mg_dm_total/sum(mg_dm_total)) %>% 
  arrange(-number_prop) %>% 
  left_join(n_fish) 



proportions_by_fishspecies %>% 
  mutate(order = mg_prop) %>% 
  pivot_longer(cols = ends_with("prop")) %>% 
  filter(name != "n_prop") %>%
  mutate(name = case_when(name == "mg_prop" ~ "% biomass",
                          name == "n_prop" ~ "Occurrence",
                          TRUE ~ "% abundance"),
         name = fct_relevel(name, "% biomass")) %>% 
  filter(value != 0) %>% 
  ggplot(aes(x = reorder(prey_taxon, order), y = value*100, color = name)) +
  geom_point() +
  scale_color_colorblind() +
  geom_line(aes(group = name)) +
  coord_flip() +
  labs(x = "Prey Taxa (ranked by biomass)",
       y = "% in fish diets") +
  scale_y_log10() +
  theme_bw() +
  theme(legend.title = element_blank()) +
  facet_wrap(~fish_species) +
  NULL




proportions_by_fishspecies %>% 
  mutate(order = mg_prop) %>% 
  pivot_longer(cols = ends_with("prop")) %>% 
  filter(name != "n_prop") %>%
  mutate(name = case_when(name == "mg_prop" ~ "% biomass",
                          name == "n_prop" ~ "Occurrence",
                          TRUE ~ "% abundance"),
         name = fct_relevel(name, "% biomass")) %>% 
  filter(value != 0, 
         name != "% abundance") %>% 
  left_join(fish_guilds) %>% 
  ggplot(aes(x = reorder(prey_taxon, order), y = value*100)) +
  # geom_point() +
  geom_bar(stat = "identity", aes(group = fish_species, fill = prey_taxon), position = "dodge") +
  # scale_color_colorblind() +
  # geom_line(aes(group = fish_species)) +
  coord_flip() +
  labs(x = "Prey Taxa (ranked by biomass)",
       y = "% in fish diets") +
  # scale_y_log10() +
  theme_bw() +
  theme(legend.title = element_blank()) +
  facet_wrap(~fish_guild) +
  NULL


taxa_abund <- proportions_by_fishspecies %>% 
  mutate(order = mg_prop) %>% 
  pivot_longer(cols = ends_with("prop")) %>% 
  filter(name != "n_prop") %>%
  mutate(name = case_when(name == "mg_prop" ~ "% biomass",
                          name == "n_prop" ~ "Occurrence",
                          TRUE ~ "% abundance"),
         name = fct_relevel(name, "% biomass")) %>% 
  filter(value != 0, 
         name != "% biomass") %>% 
  left_join(fish_guilds) %>% 
  ggplot(aes(x = reorder(prey_taxon, order), y = value*100)) +
  # geom_point() +
  geom_bar(stat = "identity", 
           aes(group = fish_species, fill = prey_taxon), position = "dodge") +
  # scale_color_colorblind() +
  # geom_line(aes(group = fish_species)) +
  coord_flip() +
  labs(x = "Prey Taxa (ranked by biomass)",
       y = "% in fish diets",
       title = "a) % abundance") +
  # scale_y_log10() +
  theme_bw() +
  theme(legend.title = element_blank()) +
  facet_wrap(~fish_guild) +
  guides(fill = F) +
  NULL


taxa_mass <- proportions_by_fishspecies %>% 
  mutate(order = mg_prop) %>% 
  pivot_longer(cols = ends_with("prop")) %>% 
  filter(name != "n_prop") %>%
  mutate(name = case_when(name == "mg_prop" ~ "% biomass",
                          name == "n_prop" ~ "Occurrence",
                          TRUE ~ "% abundance"),
         name = fct_relevel(name, "% biomass")) %>% 
  filter(value != 0, 
         name != "% abundance") %>% 
  left_join(fish_guilds) %>% 
  ggplot(aes(x = reorder(prey_taxon, order), y = value*100)) +
  # geom_point() +
  geom_bar(stat = "identity", 
           aes(group = fish_species, fill = prey_taxon), position = "dodge") +
  # geom_line(aes(group = fish_species)) +
  coord_flip() +
  labs(x = "Prey Taxa (ranked by biomass)",
       y = "% in fish diets",
       title = "b) % dry mass") +
  # scale_y_log10() +
  theme_bw() +
  theme(legend.title = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()) +
  facet_wrap(~fish_guild) +
  guides(fill = F) +
  NULL


library(patchwork)

taxa_plot <- taxa_abund + taxa_mass

ggsave(taxa_plot, file = "plots/taxa_plot.jpg", width = 8, height = 8, dpi = 500)



proportions_by_fishspecies %>% 
  mutate(order = mg_prop) %>% 
  pivot_longer(cols = ends_with("prop")) %>% 
  filter(name != "n_prop") %>%
  mutate(name = case_when(name == "mg_prop" ~ "% biomass",
                          name == "n_prop" ~ "Occurrence",
                          TRUE ~ "% abundance"),
         name = fct_relevel(name, "% biomass"),
         chiro_other = case_when(prey_taxon == "Chironomidae" ~ "Chironomidae",
                                 TRUE ~ "Other Taxa")) %>% 
  ggplot(aes(x = reorder(fish_species, -value), y = value*100, color = chiro_other)) +
  geom_point() +
  scale_color_colorblind() +
  geom_line(aes(group = name)) +
  coord_flip() +
  labs(x = "Prey Taxa (ranked by biomass)",
       y = "% in fish diets") +
  # scale_y_log10() +
  theme_bw() +
  theme(legend.title = element_blank()) +
  NULL


# which prey was max?
proportions_by_fishspecies %>% 
  group_by(fish_species) %>% 
  filter(mg_prop == max(mg_prop)) %>% 
  left_join(fish_guilds) %>% View()


proportions_by_fishspecies %>% 
  left_join(fish_guilds) %>% 
  filter(prey_taxon == "Chironomidae") %>% 
  ggplot(aes(x = fish_guild, y = mg_prop)) +
  geom_boxplot(aes(group = fish_guild), width = 0.3) + 
  geom_point(position = position_jitter(width = 0.1)) 


#proportions by prey_class
guild_diet_multi_drymass %>% 
  select(prey_class, prey_family, prey_taxon, prey_stage, prey_ecosystem, number, sample_mg_dm) %>%
  group_by(prey_class) %>% 
  summarize(number_total = sum(number),
            mg_dm_total = sum(sample_mg_dm)) %>% 
  ungroup() %>% 
  mutate(number_prop = number_total/sum(number_total),
         mg_prop = mg_dm_total/sum(mg_dm_total)) %>% 
  arrange(-number_prop) %>% 
  left_join(n_class) %>% View()


#check for relationship between fish length and mass in sample
total_mass_permm <- guild_diet_multi_drymass %>% 
  group_by(date, fish_guild, site, fish_species, sample_id, length_mm) %>% 
  summarize(sample_mg_dm01 = sum(sample_mg_dm) + 0.01) %>% 
  ungroup() %>% 
  mutate(length_mm_number = parse_number(length_mm),
         sample_mg_dm01_permm = sample_mg_dm01/length_mm_number)

sample_mass_permm <- total_mass_permm %>% 
  ggplot(aes(x = parse_number(length_mm), y = sample_mg_dm01)) +
  geom_point() +
  facet_wrap(~fish_species) +
  labs(y = "Mass of diet sample (mg DM)",
       x = "Fish standard length (mm)")

ggsave(sample_mass_permm, file = "plots/sample_mass_permm.jpg", dpi = 300, width = 8, height = 8)


