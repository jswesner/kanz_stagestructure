source("code/functions.r") # functions
source("code/make_data.r") # make data and load packages

# load posteriors 
prop_posts <- readRDS("posteriors/prop_posts.rds")

# data to plot
prop_data <- readRDS(file = "data/prop_data.rds")


# Plot Species Means --------------------------------------------------------------

fish_species_averages <- prop_posts %>% 
  ggplot(aes(x = reorder(genus_species, -median), y = value)) + 
  geom_boxplot(aes(group = interaction(fish_species, name2), fill = name2), outlier.shape = NA) +
  scale_fill_grey(end = 1, start = 0.3) +
  facet_wrap(~name2, ncol = 1) +
  guides(fill = "none") +
  theme_default() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4),
        strip.text.x = element_text(hjust = -0.01)) + 
  labs(x = "Fish Species or Family",
       y = "Proportion of diet") +
  geom_jitter(data = prop_data,
              width = 0.1, height = 0, size = 0.4,
              shape = 21) +
  NULL

saveRDS(fish_species_averages, file = "plots/fish_species_averages.rds")
ggsave(fish_species_averages, file = "plots/fish_species_averages.jpg", dpi = 500, width = 6.5, height = 8.5)
ggsave(fish_species_averages, file = "plots/Fig2.eps", dpi = 600, width = 6.5, height = 8.5)

# Summarize ---------------------------------------------------------------
library(tidybayes)
library(scales)

# proportion nonconsumers
prop_posts %>% 
  filter(fish_species == "Average") %>% 
  group_by(name) %>% 
  median_qi(value)

prop_posts %>% 
  filter(fish_species != "Average") %>% 
  group_by(name, fish_species) %>% 
  median_qi(value) %>% 
  print(n = Inf)
  
table_summary <- prop_posts %>% 
  filter(fish_species != "Average") %>% 
  group_by(name, fish_species) %>% 
  summarize(median = round(median(value), 2),
            sd = round(sd(value), 2)) %>% 
  mutate(mean_sd = paste0(median, " (", sd, ")")) %>% 
  select(-median, -sd) %>% 
  pivot_wider(names_from = name, values_from = mean_sd) %>% 
  arrange(desc(prop_noncons)) %>% 
  select(fish_species, prop_noncons, prop_noncons_aquatic, prop_terrestrial, prop_chiro_pup_adult)
  
write_csv(table_summary, file = "tables/table_summary.csv")  
