######################################################################################
## preprocess_experimental_data.R
##
## Calculates the summary statistics of the experimental data
##
## Author: Jana S. Huisman
## Last update: Feb 2021
######################################################################################
source("functions/summary_statistics.R")

library(tidyverse)
library("ggplot2")

################

data_folder = "../data/"
plot_folder = "../figures/"

################
# Read in the experimental data ####

mouse_subset = 1:13

# all in feces
pop_timecourse_data_raw = read_csv(paste0(data_folder, 'pop_sizes_feces.csv'))
pop_timecourse_data <- pop_timecourse_data_raw %>%
  pivot_longer(cols = 2:14, names_to = 'Mouse', values_to = 'Pop_size') %>%
  filter(Mouse %in% mouse_subset)

all_tags_raw = read_csv(paste0(data_folder, "all_tags_feces.csv"))

# each row sums to one
all_tags_long = all_tags_raw %>%
  rowwise() %>%
  mutate(chromo_tot = sum(chromo_tag_1, chromo_tag_2, chromo_tag_3),
         plasmid_tot = sum(plasmid_tag_1, plasmid_tag_2, plasmid_tag_3, plasmid_tag_4)) %>%
  pivot_longer(cols = 2:8, names_to = 'Tag') %>%
  mutate(Type = ifelse(grepl('chromo', Tag), 'Chromo', 'Plasmid')) %>%
  mutate(norm_value = ifelse(Type == 'Chromo', value/chromo_tot, value/plasmid_tot)) %>%
  filter(Mouse %in% mouse_subset)

###########################################################
# compute summary statistics ####

reseeding_timing <- pop_timecourse_data %>%
  mutate(reseeding = Pop_size >= 1e6) %>%
  group_by(Type, Mouse) %>% 
  summarise(first_day = ifelse(any(reseeding), Day[min(which(reseeding))], 11),
            .groups = 'drop') %>%
  mutate(first_day = ifelse(is.na(first_day), 11, first_day)) %>%
  group_by(Type) %>%
  summarise(mean = mean(first_day),
            sd = sd(first_day),
            .groups = 'drop')

final_pop_size <- pop_timecourse_data %>%
  filter(Day == 10) %>%
  group_by(Type) %>%
  summarise(mean = mean(Pop_size),
            sd = sd(Pop_size),
            .groups = 'drop')

tag_dist <- all_tags_long %>%
  filter(!is.na(value)) %>%
  group_by(Type, Mouse) %>%
  summarise(ratio = tag_ratio(norm_value, detection_limit=8.9e-5),
            evenness = evenness(norm_value, detection_limit = 8.9e-5),
            prop = sum(norm_value > 8.9e-5)/n(),
            .groups = 'drop') %>%
  pivot_longer(cols = c(ratio, evenness, prop), names_to = "measure") %>%
  group_by(Type, measure) %>%
  summarise(mean = mean(value, na.rm = T),
            sd = sd(value, na.rm = T),
            .groups = 'drop')

# save experimental summary statistics ####

exp_stats <- bind_rows(reseeding_timing %>%
              mutate(measure = 'Reseeding timing',
              mean = mean - 2), # to make it comparable to the simulation results
          final_pop_size %>%
            mutate(measure = 'Final Pop Size'),
          tag_dist)

write_csv(exp_stats, paste0(data_folder, 'exp_stats.csv'))




