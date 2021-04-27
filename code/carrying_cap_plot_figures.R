######################################################################################
## plot_carrying_cap_figures.R
##
##
## Author: Jana S. Huisman
## Last update: Feb 2021
######################################################################################
# change the commented lines depending on whether you're plotting results for
# optimal, higher, or lower conjugation rates
######################################################################################

library("ggplot2")
library("tidyverse")
library(viridis)

################

#simulation_data_folder = "../simulations_cc/"
#simulation_data_folder = "../simulations_cc_high_conj/"
simulation_data_folder = "../simulations_cc_low_conj/"
plot_folder = "../figures/"

################
ensemble_size = 100
carrying_cap_options = 10**seq(3,10,0.25)

####### Load Sim Data #########

all_simulations <- data.frame()
for (carrying_cap_i in 1:length(carrying_cap_options)){
    sim_data <- try(read_csv(paste0(simulation_data_folder, carrying_cap_i, '.csv'),
                             col_types = cols(
                               X1 = col_double(),
                               time_CFU_trans = col_double(),
                               time_CFU_recip = col_double(),
                               Even_trans = col_double(),
                               prop_trans = col_double(),
                               ratio_trans = col_double(),
                               Even_recip = col_double(),
                               prop_recip = col_double(),
                               ratio_recip = col_double(),
                               trans_pop_size = col_double(),
                               recip_pop_size = col_double()
                             )) )
    if ('try-error' %in% class(sim_data)){
      next
    }
    
    all_simulations <- sim_data %>%
      rename(replicate = X1) %>%
      mutate(carrying_cap = carrying_cap_options[carrying_cap_i]) %>%
      bind_rows(all_simulations)
    
}

#write_csv(all_simulations, paste0(simulation_data_folder, 'all_simulation_results.csv'))

####### Plot Sim Data #########

plot_df <- all_simulations %>%
  select(carrying_cap, replicate, Even_recip, Even_trans, recip_pop_size, trans_pop_size,
         time_CFU_recip, time_CFU_trans) %>%
  mutate(recip_pop_size = log10(recip_pop_size), 
         trans_pop_size = log10(trans_pop_size)) %>%
  pivot_longer(cols = 3:8, names_to = 'measure') %>%
  mutate(measure = factor(measure, levels = c('Even_recip', 'Even_trans', 'recip_pop_size', 'trans_pop_size', 
                                              'time_CFU_recip', 'time_CFU_trans'))) %>%
  filter(measure %in% c('recip_pop_size', 'trans_pop_size') )

method_names <- setNames(
  c("Evenness Chromosomal Tags",
    "Evenness Plasmid Tags",
    "Recipient Pop. Size \n(log10, day 10)",
    "Transconjugant Pop. Size \n(log10, day 10)",
    "Time to 10^6 CFU/g feces \nRecipients",
    "Time to 10^6 CFU/g feces \nTransconjugants"),
  c('Even_recip', 'Even_trans', 'recip_pop_size', 'trans_pop_size', 
    'time_CFU_recip', 'time_CFU_trans'))

ggplot(plot_df) +
  geom_boxplot(aes(x = log10(carrying_cap), 
                y = value, group = log10(carrying_cap) )) +
  facet_wrap(vars(measure), ncol =2, scale = 'free_y', 
             labeller = labeller(measure = method_names)) +
  scale_colour_viridis() + 
  labs( x = 'Carrying Capacity (log10)', 
        y = 'Measure') +
  theme_minimal() +
  theme(
    strip.text.x= element_text(size=25),
    strip.text.y= element_text(size=25),
    axis.text.y= element_text(size=20),
    axis.text.x= element_text(size=20),
    axis.title.y =  element_text(size=25),
    axis.title.x =  element_text(size=25),
    legend.text= element_text(size=20),
    legend.title= element_text(size=25)
  )

#ggsave(paste0(plot_folder, 'carrying_cap.pdf'), width = 12, height = 15)
#ggsave(paste0(plot_folder, 'carrying_cap_high_conj.pdf'), width = 12, height = 15)
ggsave(paste0(plot_folder, 'carrying_cap_low_conj.pdf'), width = 12, height = 15)

df_to_save <- all_simulations %>%
  select(carrying_cap, replicate, recip_pop_size, trans_pop_size)
#write_csv(df_to_save, paste0(simulation_data_folder, "carrying_cap.csv"))
#write_csv(df_to_save, paste0(simulation_data_folder, "carrying_cap_high_conj.csv"))
write_csv(df_to_save, paste0(simulation_data_folder, "carrying_cap_low_conj.csv"))

