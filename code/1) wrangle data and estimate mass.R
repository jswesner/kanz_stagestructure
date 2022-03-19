library(tidyverse)
library(readxl)
library(googledrive)
library(janitor)
library(EBImage)
library(pdftools)
library(taxize)
library(ggridges)
library(readxl)
library(lubridate)



##################################################################################################
#########   This script converts raw data to the curated dataset 'guild_diet_multi_drymass'
#########   To skip these steps (which require internet and an API from taxize), use this code

guild_diet_multi_drymass <- readRDS(file = "data/guild_diet_multi_drymass.rds")

#########   To re-create the final dataset by hand, follow the code below
##################################################################################################

#get taxized list from fish database
prey_taxa_have <- read_csv("data/data_to_compile/prey_taxa_have.csv") %>% 
  mutate(prey_taxon = str_to_sentence(prey_taxon)) %>% 
  distinct(prey_taxon, .keep_all = T) %>% 
  select(-life_stage, -superclass, -subclass, -subphylum, -type)


# load full dataset
full_data_taxized <- read_excel("data/data_to_compile/NEW_diet2019_fixed_images.xlsx") %>%
  type_convert() %>% 
  pivot_longer(cols = c(-pic_1,
                        -pic_2,
                        -pic_3,
                        -pic_4,
                        -pic_5,
                        -pic_6,
                        -sample_id,
                        -site,
                        -date,
                        -fish,
                        -entered_by,
                        -notes,
                        -id,
                        -length_mm))  %>% 
  rename(prey_type = name) %>% 
  separate(prey_type, c("prey_taxon", "prey_stage", "prey_ecosystem")) %>% 
  mutate(prey_taxon = str_to_sentence(prey_taxon),
         life_stage = case_when(prey_stage == "a" ~ "adult",
                                prey_stage == "l" ~ "larvae",
                                prey_stage == "p" ~ "pupae",
                                TRUE ~ prey_stage)) %>% 
  left_join(prey_taxa_have, by = "prey_taxon") %>% 
  mutate(source = "gut_contents")


#load length data into R
diet2019_with_lengths <- read_csv("data/data_to_compile/diet2019_with_lengths.csv") %>% 
  rename(prey_type = name)

d_taxized <- diet2019_with_lengths %>% 
  pivot_longer(cols= starts_with("mm_"), 
               names_to = "measure",
               values_to = "prey_length") %>% 
  separate(prey_type, c("prey_taxon", "prey_stage", "prey_ecosystem")) %>% 
  mutate(prey_taxon = str_to_sentence(prey_taxon),
         life_stage = case_when(prey_stage == "a" ~ "adult",
                                prey_stage == "l" ~ "larvae",
                                prey_stage == "p" ~ "pupae",
                                TRUE ~ prey_stage)) %>% 
  left_join(prey_taxa_have, by = "prey_taxon") %>% 
  mutate(source = "gut_contents")





# Length-mass regression equations ----------------------------------

length_mass_equations <- read_csv("data/data_to_compile/length_mass_regressions.csv") %>% 
  separate(b, c("b_mean", "b_se"), sep = "_") %>% 
  separate(a, c("a_mean", "a_se"), sep = "_") %>% 
  filter(mass_units == "mg dry mass") %>% 
  select(-taxon) %>% 
  type_convert() %>% 
  rename(length_measure = measure) %>% 
  unite("prey_orderstage", c(order, life_stage)) %>% 
  distinct(a_mean, b_mean, prey_orderstage)

d_taxized_length <- d_taxized %>% 
  unite("prey_orderstage", c(order, life_stage)) %>% 
  distinct(prey_orderstage, prey_length) %>% 
  left_join(length_mass_equations, by = "prey_orderstage") %>% 
  mutate(mg_dm_ind = a_mean*prey_length^b_mean ) %>% 
  separate(prey_orderstage, c("order", "life_stage")) %>% 
  mutate(source = "gut_contents")

# masses from bow creek
mass_data_bowcreek <- read_excel("data/data_to_compile/Bow Creek Insect weights.xlsx") %>% clean_names() %>%
  filter(number_of_insects >0) %>%
  drop_na(number_of_insects) %>%
  mutate(mg_ind = mg/number_of_insects,
         life_stage = case_when(grepl("arvae", life_stage) ~ "larvae",
                                grepl("upae", life_stage) ~ "pupae",
                                grepl("ymph", life_stage) ~ "larvae",
                                grepl("dult", life_stage) ~ "adult",
                                TRUE ~ life_stage)) %>%
  drop_na(mg_ind) %>%
  filter(sample == "Dry",
         !is.na(life_stage)) %>%
  mutate(source = "Warmbold 2016 - MS Thesis - USD") %>%
  rename(fam = family,
         ord = order)

# taxize
# bow_taxa <- unique(mass_data_bowcreek$lowest_taxa)

# bow_classified <- classification(bow_taxa, db = "gbif", return_id = T)

# classified <- rbind(bow_classified) %>% distinct(name, rank, query) %>%
# as_tibble() %>%
# rename(lowest_taxa = query)

# bow_taxized <- left_join(mass_data_bowcreek %>% distinct(lowest_taxa), classified) %>% 
# pivot_wider(names_from = rank, values_from = name) %>% 
# select(-`NA`)

# saveRDS(bow_taxized, file = "data/data_to_compile/bow_taxized.rds")
bow_taxized <- readRDS("data/data_to_compile/bow_taxized.rds")

#add to full dataset

mass_data_bowcreek_taxized <- left_join(mass_data_bowcreek, bow_taxized) %>% 
  select(life_stage, source, kingdom, family, phylum, class, order, genus, mg_ind) 

# compare diet estimates with larger dataset below from bow creek
d_taxized_length %>% 
  select(order, life_stage, mg_dm_ind) %>% 
  bind_rows(mass_data_bowcreek_taxized %>% select(source, order, life_stage, mg_ind) %>% rename(mg_dm_ind = mg_ind)) %>% 
  filter(!is.na(mg_dm_ind)) %>% 
  ggplot(aes(x = order, y = mg_dm_ind, color = source)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.5,
             jitter.width = 0.1), alpha = 0.2) +
  geom_boxplot(aes(group = interaction(source,order))) +
  facet_wrap(~life_stage, scales = "free") +
  scale_y_log10() +
  coord_flip()




all_mass <- bind_rows(mass_data_bowcreek_taxized, d_taxized_length %>% rename(mg_ind = mg_dm_ind) %>% 
                        filter(life_stage != "pupae"))  

# by order and life stage
order_life_stage_meansd <- all_mass %>% 
  group_by(order, life_stage) %>% 
  summarize(mean_dm = mean(mg_ind, na.rm = T),
            sd_dm = sd(mg_ind, na.rm = T))



temp_data_with_mass <- full_data_taxized %>% 
  left_join(order_life_stage_meansd)

have <- temp_data_with_mass %>% filter(!is.na(mean_dm))
still_need <- temp_data_with_mass %>% filter(is.na(mean_dm))

# by_class_only
class_meansd <- all_mass %>% 
  group_by(class) %>% 
  summarize(mean_dm = mean(mg_ind, na.rm = T),
            sd_dm = sd(mg_ind, na.rm = T))

add_class <- still_need %>% 
  select(-mean_dm, -sd_dm) %>% 
  left_join(class_meansd)

have2 <- have %>% bind_rows(add_class) 

still_need2 <- have2 %>% filter(is.na(mean_dm))

# by_phylum_only
phylum_meansd <- all_mass %>% 
  group_by(phylum) %>% 
  summarize(mean_dm = mean(mg_ind, na.rm = T),
            sd_dm = sd(mg_ind, na.rm = T))

add_phylum <- still_need2 %>% 
  select(-mean_dm, -sd_dm) %>% 
  left_join(phylum_meansd)

have3 <- have2 %>% filter(!is.na(mean_dm)) %>% 
  bind_rows(add_phylum)

still_need3 <- have3 %>% filter(is.na(mean_dm))


still_need3 %>% distinct(prey_taxon)


# load feeding categories for prey

feeding_categories <- read_csv("data/data_to_compile/feeding_categories.csv") %>% select(prey_taxon, prey_feeding, prey_stage)

# Make final diet data set ------------------------------------------------


#load fish data

fish_taxa <- readRDS("data/data_to_compile/guild_diet_multi.rds") %>% distinct(fish, species, finalguild) %>% 
  rename(fish_species = species,
         fish_guild = finalguild)

#make data set
guild_diet_multi_drymass <- have3 %>% 
  mutate(mean_dm = case_when(prey_taxon == "Crayfish" ~ 80,
                             prey_taxon == "Frog" ~ mean(c(50, 56)),
                             prey_taxon == "Branchiopoda" ~ 0.28,
                             prey_taxon == "Chironomidae" & prey_stage == "p" ~ 0.58,
                             TRUE ~ mean_dm),
         sd_dm = case_when(prey_taxon == "Frog" ~ sd(c(50,56)),
                           prey_taxon == "Branchiopoda" ~ 0.68,
                           prey_taxon == "Chironomidae" & prey_stage == "p" ~ 0.44,
                           TRUE ~ sd_dm),
         source = case_when(prey_taxon== "Frog" ~ "Kanz and Wesner - direct measure of dry mg",
                            prey_taxon == "Crayfish" ~ "Kanz and Wesner - direct measure of dry mg",
                            prey_taxon == "Branchiopoda" ~ "Robson et al. 2016 shiny app",
                            prey_taxon == "Chironomidae" & prey_stage == "p" ~ "Seidel and Wesner unpublished + 1 measure from this study",
                            TRUE ~ source),
         class = case_when(prey_taxon == "Caelifera" ~ "Insecta", TRUE ~ class),
         order = case_when(prey_taxon == "Caelifera" ~ "Orthoptera", TRUE ~ order),
         order = case_when(prey_taxon == "Empididae" ~ "Diptera", TRUE ~ order),
         prey_ecosystem = case_when(family == "Hydropsychidae" ~ "aquatic",
                                    family == "Haliplidae" ~ "aquatic",
                                    family == "Gyrinidae" ~ "aquatic",
                                    prey_taxon == "Frog" ~ "aquatic",
                                    TRUE ~ prey_ecosystem)) %>% 
  rename(prey_kingdom = kingdom,
         prey_phylum = phylum,
         prey_class = class,
         prey_order = order,
         prey_family = family,
         prey_species = species,
         number = value,
         mean_individual_drymass = mean_dm,
         sd_individual_drymass = sd_dm) %>% 
  select(-taxon) %>% 
  mutate(date = as.factor(date),
                       number = as.integer(number),
                       prey_taxaclassstage = paste0(prey_class,"_",prey_taxon,"_", prey_stage)) %>% 
  mutate(sample_mg_dm = number*mean_individual_drymass) %>% 
  left_join(fish_taxa) %>% 
  left_join(feeding_categories) %>%
  mutate(number = replace_na(number, 0),
         sample_mg_dm = replace_na(sample_mg_dm, 0)) %>% 
  mutate(prey_feeding_overall = prey_feeding,
         prey_feeding = case_when(prey_ecosystem == "terrestrial" ~ "non_consumer",
                                  is.na(prey_feeding) ~ "consumer", 
                                  prey_taxon %in% c("Dolichopodidae", "Culicidae", "Empididae", "Simuliidae",
                                                    "Zygoptera", "Anisoptera") & prey_stage == "a" ~ "non_consumer",
                                  
                                  TRUE ~ prey_feeding)) %>% 
  mutate(sample_mg_dm01 = sample_mg_dm + 0.01,
         sample_mg_dm01_permm = sample_mg_dm01/parse_number(length_mm)) %>% 
  mutate(prey_ecosystem = case_when(is.na(prey_ecosystem) ~ "unknown",
                                    TRUE ~ prey_ecosystem)) %>% 
  mutate(sample_mg_dm01 = sample_mg_dm + 0.01,
         sample_mg_dm01_permm = sample_mg_dm01/parse_number(length_mm)) %>% 
  mutate(prey_ecosystem = case_when(is.na(prey_ecosystem) ~ "unknown",
                                    TRUE ~ prey_ecosystem),
         date_formatted = ymd(date),
         date_no = yday(date_formatted) - min(yday(date_formatted)) + 1) 

chiros_only <- guild_diet_multi_drymass %>% 
  filter(grepl("hironomid", prey_family)) %>% 
  ungroup() %>% 
  mutate(mean_date_no = mean(date_no),
         original_date_no = date_no,
         date_no = date_no - mean(date_no)) %>%
  mutate(fish_prey_stage = paste0(fish_species,"_", prey_stage)) %>% 
  select(date, mean_date_no, original_date_no, date_no, 
         sample_id, fish_prey_stage, sample_mg_dm, sample_mg_dm01, site, fish_species, prey_stage)

write_csv(chiros_only, file = "data/chiros_only.csv")
saveRDS(chiros_only, file = "data/chiros_only.rds")
saveRDS(guild_diet_multi_drymass, file = "data/guild_diet_multi_drymass.rds")
write_csv(guild_diet_multi_drymass, file = "data/guild_diet_multi_drymass.csv")





