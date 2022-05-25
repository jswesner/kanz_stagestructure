
compare_proportions_df <- prey_terr_noncons %>% 
  ungroup() %>% 
  group_by(prey_feeding, site, date, fish_species, sample_id) %>% 
  summarize(value_mg = sum(sample_mg_dm, na.rm = T),
            value_no = sum(sample_number, na.rm = T)) %>% 
  group_by(sample_id) %>% 
  mutate(total_mg = sum(value_mg, na.rm = T),
         total_no = sum(value_no, na.rm = T),
         proportion_mg = value_mg/total_mg,
         proportion_no = value_no/total_no) %>% 
  ungroup() %>% 
  arrange(proportion_mg) %>% 
  mutate(rank = 1:nrow(.)) %>% 
  group_by(fish_species, prey_feeding) %>% 
  mutate(median = median(proportion_mg, na.rm = T)) %>% 
  pivot_longer(cols = c(proportion_mg, proportion_no)) 
  
library(viridis)
compare_proportions_df %>% 
  filter(prey_feeding == "non_consumer") %>% 
  ggplot(aes(y = reorder(fish_species, median), x = value, fill = median)) + 
  geom_density_ridges(aes(group = interaction(fish_species, name)),
                      alpha = 0.7,
                      point_shape = "|",
                      points = T) +
  facet_wrap(~name) +
  scale_fill_viridis() +
  coord_cartesian(xlim = c(0, 1)) +
  NULL

compare_proportions_df %>% 
  ggplot(aes(x = rank, y = value, color = name)) + 
  geom_point()

