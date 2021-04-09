######################################################################################
## plot_simulation_figures.R
##
## Calculates the summary statistics of the simulations
##
## Author: Jana S. Huisman
## Last update: Feb 2021
######################################################################################

library("ggplot2")
library("tidyverse")
library(viridis)

################

simulation_data_folder = "../simulations/"
data_folder = "../data/"
plot_folder = "../figures/"

################
ensemble_size = 100
conj_rate_options = 10**seq(-12,-3,0.25)
mig_rate_options = 10**seq(-6,3,0.25) 

####### Load Sim Data #########

all_simulations <- data.frame()
for (conj_i in 1:length(conj_rate_options)){
  for (mig_i in 1:length(mig_rate_options)){
    sim_data <- try(read_csv(paste0(simulation_data_folder, conj_i, '_',
                    mig_i, '.csv'),
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
      mutate(conj_rate = conj_rate_options[conj_i],
             mig_rate = mig_rate_options[mig_i]) %>%
      bind_rows(all_simulations)
    
  }
}

write_csv(all_simulations, paste0(simulation_data_folder, 'all_simulation_results.csv'))

####### Plot Simulated Data Means #########

sim_means <- all_simulations %>%
  group_by(conj_rate, mig_rate) %>%
  summarise(across(everything(), ~ mean(., na.rm = T)),
            .groups = 'drop') %>%
  select(-replicate)

sim_sds <- all_simulations %>%
  group_by(conj_rate, mig_rate) %>%
  summarise(across(everything(), ~ sd(., na.rm = T)),
            .groups = 'drop') %>%
  select(-replicate)
  

method_names <- setNames(
         c("Time to \n10^6 CFU/g \nfeces \nTransconjugants",
           "Time to \n10^6 CFU/g \nfeces \nRecipients",
           "Evenness \nPlasmid \nTags",
           "Proportion \nPlasmid \nTags",
           "Ratio \n1st/2nd \nPlasmid \nTag",
           "Evenness \nChromo \nTags",
           "Proportion \nChromo \nTags",
           "Ratio \n1st/2nd \nChromo \nTag",
           "Final \nTransconjugant \nPop. Size",
           "Final \nRecipient \nPop. Size"),
         colnames(sim_means[3:12]))

## Not shown in paper ####
# for (measure in names(method_names)){
#   ggplot(sim_means) +
#     geom_tile(aes(x = log10(conj_rate), y = log10(mig_rate), 
#                   fill = get(measure) )) +
#     scale_fill_viridis() + 
#     labs( x = 'Conjugation rate (log10)', y = 'Migration rate (log10)', 
#           fill = method_names[measure]) +
#     theme_minimal() +
#     theme(
#       axis.text.y= element_text(size=20),
#       axis.text.x= element_text(size=20),
#       axis.title.y =  element_text(size=25),
#       axis.title.x =  element_text(size=25),
#       legend.text= element_text(size=20),
#       legend.title= element_text(size=25)
#     )
#   
#   ggsave(paste0(plot_folder, measure, '.pdf'))
#   
#   df_to_save <- sim_means %>%
#     select(conj_rate, mig_rate, value = measure) %>%
#     pivot_wider(names_from = conj_rate, values_from = value)
#   write_csv(df_to_save,
#             paste0(simulation_data_folder, measure,".csv"))
#   
#   ## Variance
#   ggplot(sim_sds) +
#     geom_tile(aes(x = log10(conj_rate), y = log10(mig_rate), 
#                   fill = get(measure) )) +
#     scale_fill_viridis() + 
#     labs( x = 'Conjugation rate (log10)', y = 'Migration rate (log10)', 
#           fill = method_names[measure]) +
#     theme_minimal() +
#     theme(
#       axis.text.y= element_text(size=20),
#       axis.text.x= element_text(size=20),
#       axis.title.y =  element_text(size=25),
#       axis.title.x =  element_text(size=25),
#       legend.text= element_text(size=20),
#       legend.title= element_text(size=25)
#     )
#   
#   ggsave(paste0(plot_folder, 'sd_', measure, '.pdf'))
# }

###########################################################
# Plot the values wrt. experimental means ####

exp_stats <- read_csv(paste0(data_folder, 'exp_stats.csv'))

exp_stats['summary_stat'] = c("time_CFU_donor", "time_CFU_recip", "time_CFU_trans",
                              "donor_pop_size", "recip_pop_size", "trans_pop_size", 
                              "Even_recip",     "prop_recip",     "ratio_recip",
                              "Even_trans",     "prop_trans",     "ratio_trans")

for (measure in names(method_names)){
  ggplot(sim_means) +
    geom_tile(aes(x = log10(conj_rate), y = log10(mig_rate), 
                  fill = get(measure) )) +
    scale_fill_gradient2(high = 'red', low = 'blue', midpoint =  exp_stats[exp_stats$summary_stat == measure, ]$mean ) +
    labs( x = 'Conjugation rate (log10)', y = 'Migration rate (log10)', 
          fill = method_names[measure]) +
    theme_minimal() +
    theme(
      axis.text.y= element_text(size=20),
      axis.text.x= element_text(size=20),
      axis.title.y =  element_text(size=25),
      axis.title.x =  element_text(size=25),
      legend.text= element_text(size=20),
      legend.title= element_text(size=25)
    )
  
  ggsave(paste0(plot_folder, measure, '_exp.pdf'))

}


