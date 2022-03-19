source("code/make_data.r")

library(RInSp)
library(tidybayes)

# load posts just to sort by median propnoncosumers later
posts_cons_noncons_fishspecies <- readRDS("posteriors/posts_all.rds") %>% 
  group_by(.draw, fish_species, prey_feeding) %>% 
  summarize(total = sum(.epred)) %>% 
  pivot_wider(names_from = prey_feeding, values_from = total) %>% 
  mutate(total = consumer + non_consumer,
         prop_nonconsumer = non_consumer/total) %>% 
  group_by(fish_species)  %>% 
  median_qi() %>% 
  arrange(fish_species)

# raw diet data without prey stages
prey_taxa_diet <- guild_diet_multi_drymass %>% 
  group_by(sample_id, site, date, fish_species) %>% 
  mutate(total = sum(sample_mg_dm)) %>% 
  filter(total > 0) %>% 
  select(sample_id, site, date, fish_species, sample_mg_dm, prey_taxon, prey_feeding) %>% 
  group_by(sample_id, site, date, fish_species, prey_taxon) %>% 
  summarize(sample_mg_dm = sum(sample_mg_dm)) %>% 
  pivot_wider(names_from = prey_taxon, values_from = sample_mg_dm)

# Use RInSp to calculate overlap for individuals from Bolnick et al. (2002)
burbank_rinsp <- import.RInSp(prey_taxa_diet %>% filter(site == "burbank"), row.names = 1, info.cols = c(1:4))
gunderson_rinsp <- import.RInSp(prey_taxa_diet %>% filter(site == "gunderson"), row.names = 1, info.cols = c(1:4))
spiritmound_rinsp <- import.RInSp(prey_taxa_diet %>% filter(site == "spirit mound"), row.names = 1, info.cols = c(1:4))
littlebridge_rinsp <- import.RInSp(prey_taxa_diet %>% filter(site == "little bridge"), row.names = 1, info.cols = c(1:4))

# extract overlap
burbank_similarity <- overlap(burbank_rinsp)
gunderson_similarity <- overlap(gunderson_rinsp)
spiritmound_similarity <- overlap(spiritmound_rinsp)
littlebridge_similarity <- overlap(littlebridge_rinsp)

# convert to tibble
burbank_similarity <- as_tibble(burbank_similarity)
gunderson_similarity <- as_tibble(gunderson_similarity)
spiritmound_similarity <- as_tibble(spiritmound_similarity)
littlebridge_similarity <- as_tibble(littlebridge_similarity)

# append overlap values to individual samples
burbank_overlap_taxafeeding <- burbank_rinsp$info %>% as_tibble() %>% 
  mutate(overlap = burbank_similarity$meanindividualoverlap[,1]) %>% clean_names()
gunderson_overlap_taxafeeding <- gunderson_rinsp$info %>% as_tibble() %>% 
  mutate(overlap = gunderson_similarity$meanindividualoverlap[,1]) %>% clean_names()
spiritmound_overlap_taxafeeding <- spiritmound_rinsp$info %>% as_tibble() %>% 
  mutate(overlap = spiritmound_similarity$meanindividualoverlap[,1]) %>% clean_names()
littlebridge_overlap_taxafeeding <- littlebridge_rinsp$info %>% as_tibble() %>% 
  mutate(overlap = littlebridge_similarity$meanindividualoverlap[,1]) %>% clean_names()

# combine into one data frame
all_taxa_overlap <- bind_rows(burbank_overlap_taxafeeding,
                              gunderson_overlap_taxafeeding,
                              spiritmound_overlap_taxafeeding,
                              littlebridge_overlap_taxafeeding) %>% 
  mutate(aggregation = "prey_taxa_only")



# repeat above steps but now for prey with stage assigned
prey_taxastage_diet <- guild_diet_multi_drymass %>% 
  select(sample_id, site, date, fish_species, sample_mg_dm, prey_taxon, prey_stage) %>% 
  group_by(sample_id, site, date, fish_species, prey_taxon, prey_stage) %>% 
  summarize(sample_mg_dm = sum(sample_mg_dm)) %>% 
  unite("prey_taxa_stage", c(prey_taxon, prey_stage), sep = "_") %>% 
  pivot_wider(names_from = prey_taxa_stage, values_from = sample_mg_dm)

burbankstage_rinsp <- import.RInSp(prey_taxastage_diet %>% filter(site == "burbank"), row.names = 1, info.cols = c(1:4))
gundersonstage_rinsp <- import.RInSp(prey_taxastage_diet %>% filter(site == "gunderson"), row.names = 1, info.cols = c(1:4))
spiritmoundstage_rinsp <- import.RInSp(prey_taxastage_diet %>% filter(site == "spirit mound"), row.names = 1, info.cols = c(1:4))
littlebridgestage_rinsp <- import.RInSp(prey_taxastage_diet %>% filter(site == "little bridge"), row.names = 1, info.cols = c(1:4))


burbank_similarity_stage <- overlap(burbankstage_rinsp)
gunderson_similarity_stage <- overlap(gundersonstage_rinsp)
spiritmound_similarity_stage <- overlap(spiritmoundstage_rinsp)
littlebridge_similarity_stage <- overlap(littlebridgestage_rinsp)

burbank_overlap_taxastage <- burbankstage_rinsp$info %>% as_tibble() %>% 
  mutate(overlap = burbank_similarity_stage$meanindividualoverlap[,1]) %>% clean_names()
gunderson_overlap_taxastage <- gundersonstage_rinsp$info %>% as_tibble() %>% 
  mutate(overlap = gunderson_similarity_stage$meanindividualoverlap[,1]) %>% clean_names()
spiritmound_overlap_taxastage <- spiritmoundstage_rinsp$info %>% as_tibble() %>% 
  mutate(overlap = spiritmound_similarity_stage$meanindividualoverlap[,1]) %>% clean_names()
littlebridge_overlap_taxastage <- littlebridgestage_rinsp$info %>% as_tibble() %>% 
  mutate(overlap = littlebridge_similarity_stage$meanindividualoverlap[,1]) %>% clean_names()


all_stage_overlap <- bind_rows(burbank_overlap_taxastage,
                              gunderson_overlap_taxastage,
                              spiritmound_overlap_taxastage,
                              littlebridge_overlap_taxastage) %>% 
  mutate(aggregation = "prey_taxa_stage")

# combine estimates of overlap for taxa only and taxa + stage
all_overlap <- all_stage_overlap %>% 
  bind_rows(all_taxa_overlap) 

# get raw proportion nonconsumer data for regressions later
prop_nonconsumer <- guild_diet_multi_drymass %>% 
  select(sample_id, site, date, fish_species, sample_mg_dm, prey_feeding)  %>% 
  group_by(sample_id, site, date, fish_species, prey_feeding) %>% 
  summarize(total = sum(sample_mg_dm)) %>% 
  pivot_wider(names_from = prey_feeding, values_from = total) %>% 
  mutate(total = consumer + non_consumer,
         prop_nonconsumer = non_consumer/total)

# get raw proportion non-larval data for regressions later
prop_nonlarval <- guild_diet_multi_drymass %>% 
  select(sample_id, site, date, fish_species, sample_mg_dm, prey_stage)  %>% 
  group_by(sample_id, site, date, fish_species, prey_stage) %>% 
  mutate(prey_stage = case_when(is.na(prey_stage) ~ "unknown", TRUE ~ prey_stage)) %>% 
  summarize(total = sum(sample_mg_dm)) %>% 
  pivot_wider(names_from = prey_stage, values_from = total) %>% 
  mutate(total = a + l + p + unknown,
         prop_nonlarval = (a + p)/total)

# combine results from rlinsp with prop nonconsumer/non-larval per samples
diffs <- all_overlap %>% 
  pivot_wider(names_from = aggregation, values_from = overlap) %>% 
  mutate(diff_stage_minus_taxa = prey_taxa_stage - prey_taxa_only) %>% 
  left_join(prop_nonconsumer %>% select(sample_id, prop_nonconsumer)) %>% 
  left_join(prop_nonlarval %>% select(sample_id, prop_nonlarval)) %>% 
  ungroup() %>% 
  mutate(mean_prop_nonlarval = mean(prop_nonlarval, na.rm = T),
         sd_prop_nonlarval = sd(prop_nonlarval, na.rm = T),
         prop_nonlarval_z = (prop_nonlarval - mean_prop_nonlarval)/sd_prop_nonlarval)

# make a plot of differences
diff_plot <- all_overlap %>% 
  mutate(name_yn = case_when(grepl("only", aggregation) ~ "Without life-stage", TRUE ~ "With life-stage"),
       name_yn = fct_relevel(name_yn, "Without life-stage")) %>% 
  group_by(site, aggregation) %>% 
  arrange(desc(overlap)) %>% 
  mutate(rank = row_number())  %>% 
  mutate(site = case_when(site == "burbank" ~ "a) Burbank",
                          site == "gunderson" ~ "b) Gunderson",
                          site == "little bridge" ~ "c) Little Bridge",
                          site == "spirit mound" ~ "d) Spirit Mound")) %>% 
  ggplot(aes(y = overlap, x = name_yn)) + 
  geom_point(aes(color = name_yn), size = 0.7,
             position = position_dodge(width = 0.5)) +
  geom_line(aes(group = sample_id), alpha = 0.2) +
  labs(y = "Individual Diet Overlap",
       x = "Method (with or without life-stages)",
       color = "") +
  facet_wrap(~site, ncol = 4) +
  scale_x_discrete(labels = c("With", "Without")) +
  scale_color_colorblind() +
  theme_bw() +
  guides(color = "none") +
  theme(panel.grid = element_blank(),
        legend.background = element_rect(fill = NA),
        legend.position = "top") +
  NULL

# Model diet overlap ------------------------------------------------------

# arrange data
diffs_long <- diffs %>% 
  pivot_longer(cols = c(prey_taxa_stage, prey_taxa_only)) %>% 
  select(sample_id, site, date, fish_species, value, name, prop_nonlarval, prop_nonlarval_z) %>% 
  mutate(value01 = value + 0.01)

# fit model of overlap differences
brm_dietoverlap_species <- brm(value01 ~ 1 + name*site + (1|sample_id) + 
                                    (1|date) + (1 + name*site|fish_species),
                               data = diffs_long , family = Beta(link = "logit", link_phi = "log"),
                               prior = c(prior(normal(0, 1), class = "Intercept"),
                                         prior(normal(0, 1), class = "b"),
                                         prior(exponential(2), class = "sd")),
                               file = "models/species_models/brm_dietoverlap_species.rds",
                               file_refit = "on_change", cores = 4, chains = 4, iter = 2000)

brm_dietoverlap_species <- update(brm_dietoverlap_species, iter = 2000, chains = 4)

# read model fit above
brm_dietoverlap_species <- readRDS("models/species_models/brm_dietoverlap_species.rds")

# check the model
pp_check(brm_dietoverlap_species, type = "hist")

# extract posteriors
posts_diffs_species <- brm_dietoverlap_species$data %>% 
  distinct(site, name, fish_species) %>% 
  mutate(date = NA, sample_id = NA) %>% 
  add_epred_draws(brm_dietoverlap_species, re_formula = NULL)

# summarize posteriors
posts_diffs_species %>% ungroup() %>% select(site, fish_species, name, .draw, .epred) %>% 
  # pivot_wider(names_from = name, values_from = .epred) %>% 
  group_by(site, fish_species, name) %>% 
  median_qi()

posts_diffs_species %>% ungroup() %>% select(site, fish_species, name, .draw, .epred) %>% 
  pivot_wider(names_from = name, values_from = .epred) %>%
  mutate(diff = 1 - (prey_taxa_stage/prey_taxa_only)) %>% 
  group_by(fish_species) %>%
  median_qi(diff) %>% 
  arrange(diff)

color_median <- posts_cons_noncons_fishspecies %>% 
  group_by(fish_species) %>% 
  summarize(median = median(prop_nonconsumer)) %>% 
  mutate(median = case_when(fish_species == "Average" ~ 0, T ~ median)) # places average at the bottom

# make plot of model results
species_overlap_plot <- posts_diffs_species %>% 
  left_join(color_median) %>% 
  mutate(name = fct_relevel(name, "prey_taxa_only")) %>% 
  mutate(name = case_when(name == "prey_taxa_only" ~ "With life-stages",
                          TRUE ~ "Without life-stages")) %>% 
  mutate(site = str_to_title(site)) %>% 
  ggplot(aes(x = .epred, y = reorder(fish_species, median), fill = name)) +
  # geom_density_ridges(alpha = 0.6) +
  geom_boxplot(aes(group = interaction(name, fish_species)),
  outlier.shape = NA) +
  scale_fill_grey(start = 0.6, end = 1) +
  facet_wrap(~site, nrow = 1) +
  theme_default() +
  scale_x_continuous(breaks = c(0, 0.5, 1), limits = c(0,1)) + 
  labs(fill = "Method",
       x = "Dietary Overlap") +
  theme(legend.position = "top", 
        axis.title.y = element_blank())


ggsave(species_overlap_plot, file = "plots/species_overlap_plot.jpg", width = 6, height = 6)
ggsave(species_overlap_plot, file = "plots/Fig3.eps", dpi = 600, width = 6, height = 6)


# regression of differences as function of the proportion of non-larval prey
get_prior(diff_stage_minus_taxa ~ prop_nonlarval_z + (1|date) + (1|site) + (1|fish_species),
          data = diffs , family = gaussian())

brm_dietoverlap <- brm(bf(diff_stage_minus_taxa ~ prop_nonlarval_z*site + (1|date) + (1|fish_species),
                          sigma ~ prop_nonlarval_z*site),
                       data = diffs , family = gaussian(),
                       prior = c(prior(normal(0, 1), class = "Intercept"),
                                 prior(normal(0, 1), class = "b"),
                                 prior(exponential(2), class = "sd")),
                       file = "models/species_models/brm_dietoverlap.rds",
                       file_refit = "on_change", cores = 4, chains = 4, iter = 2000)

# extract posteriors
overlap_posts <- as_draws_df(brm_dietoverlap) %>% as_tibble() %>% clean_names()

# summarize posteriors
overlap_posts %>% select(contains("b_prop_nonlarval")) %>%
  mutate(.draw = 1:nrow(.)) %>%
  pivot_longer(cols = contains("_z_")) %>% 
  mutate(group_value = b_prop_nonlarval_z + value) %>% 
  select(-value) %>% 
  pivot_wider(names_from = name, values_from = group_value) %>% 
  pivot_longer(cols = -.draw) %>% 
  mutate(site = case_when(grepl("_z_", name) ~ name,
                          TRUE ~ "burbank")) %>% 
  mutate(site = str_remove(site, "b_prop_nonlarval_z_"),
         slope = exp(value)) %>% 
  group_by(site) %>% 
  median_qi(slope)


brm_posts <- brm_dietoverlap$data %>% 
  distinct(site, prop_nonlarval_z) %>% 
  add_epred_draws(brm_dietoverlap, re_formula = NA) %>% 
  left_join(diffs %>% distinct(site, prop_nonlarval_z, prop_nonlarval))

# make plot of results
diff_regression_plot <- brm_posts %>%
  group_by(site, prop_nonlarval_z, prop_nonlarval) %>% 
  median_qi(.epred) %>% 
  mutate(site = str_to_title(site)) %>% 
  # mutate(site = case_when(site == "burbank" ~ "e) Burbank",
  #                         site == "gunderson" ~ "f) Gunderson",
  #                         site == "little bridge" ~ "g) Little Bridge",
  #                         site == "spirit mound" ~ "h) Spirit Mound")) %>% 
  ggplot(aes(x = prop_nonlarval, y = .epred, fill = site)) + 
  geom_point(data = brm_dietoverlap$data %>% 
               left_join(diffs %>% distinct(site, prop_nonlarval_z, prop_nonlarval)) %>% 
               mutate(site = str_to_title(site)),
             aes(y = diff_stage_minus_taxa, shape = site, color = site),
             size = 0.7) +
  geom_line() +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.3) +
  theme_bw() +
  scale_color_colorblind() +
  scale_fill_colorblind() +
  scale_shape_manual(values = c(1, 19, 20, 21)) +
  labs(x = "Proportion of non-larval prey in diet\n(standardized)",
       y = "Overlap differences\n('with life-stage' - 'without life-stage')",
       color = "Site", 
       fill = "Site",
       shape = "Site") +
  # facet_wrap(~site) +
  NULL


ggsave(diff_regression_plot, file = "plots/overlap_plot_single.jpg", dpi = 500, width = 6, height = 4)
ggsave(diff_regression_plot, file = "plots/Fig4.pdf", dpi = 600, width = 6, height = 4)

overlap_plot <- diff_plot/diff_regression_plot + plot_layout(heights = c(0.8, 1.5))
saveRDS(overlap_plot, file = "plots/overlap_plot.rds")
ggsave(overlap_plot, file = "plots/overlap_plot.jpg", dpi = 500, width = 6, height = 8)

overlap_plot <- readRDS("plots/overlap_plot.rds")

# summary -----------------------------------------------------------------

brm_posts %>% 
  ungroup() %>% 
  group_by(site) %>%
  median_qi(.epred)


brm_posts %>% 
  distinct(prop_nonlarval, prop_nonlarval_z) %>% 
  ggplot(aes(x = prop_nonlarval, y = prop_nonlarval_z)) + 
  geom_point()


low_mid_high <- brm_dietoverlap$data %>% 
  group_by(site) %>% 
  summarize(min = -0.01,
            median = 1,
            max = max(prop_nonlarval_z)) %>% 
  pivot_longer(cols = -site, values_to = "prop_nonlarval_z") %>% 
  add_epred_draws(brm_dietoverlap, re_formula = NA)

low_mid_high %>% 
  group_by(name) %>% 
  median_qi(.epred)

