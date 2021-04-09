######################################################################################
## compare_sim_experiment.R
##
## Calculates the summary statistics of the experimental data and compares them
## against the summary statistics of the simulations
## Plots the figures in the paper
##
## Author: Jana S. Huisman
## Last update: Feb 2021
######################################################################################
source("functions/summary_statistics.R")

library('tidyverse')
library("ggplot2")
library(viridis)

theme_set(theme_minimal() +
            theme(
              axis.text.y= element_text(size=20),
              axis.text.x= element_text(size=20),
              axis.title.y =  element_text(size=25),
              axis.title.x =  element_text(size=25),
              legend.text= element_text(size=20),
              legend.title= element_text(size=25)
            ))

################

simulation_data_folder = "../simulations/"
data_folder = "../data/"
plot_folder = "../figures/"

##################
all_simulations <- read_csv(paste0(simulation_data_folder, 'all_simulation_results.csv'),
                            col_types = 'ddddddddddddd')

exp_stats <- read_csv(paste0(data_folder, 'exp_stats.csv'))

exp_stats['summary_stat'] = c("time_CFU_donor", "time_CFU_recip", "time_CFU_trans",
                              "donor_pop_size", "recip_pop_size", "trans_pop_size", 
                              "Even_recip",     "prop_recip",     "ratio_recip",
                              "Even_trans",     "prop_trans",     "ratio_trans")


##################### Helper Functions ######################################

# Determine whether the simulation results are close enough to the experimental
# results to be accepted
acceptance_condition <- function(all_simulations, exp_stats, sd_dist = 1){
  conform_simulations <- all_simulations
  for(stat in colnames(all_simulations)[2:11]){
    target_size = exp_stats %>% filter(summary_stat == stat) %>% pull(mean)
    sd_target = exp_stats %>% filter(summary_stat == stat) %>% pull(sd)
    
    conform_simulations[stat] <- (conform_simulations[stat] < (target_size + sd_dist*sd_target) ) &
                                 (conform_simulations[stat] > (target_size - sd_dist*sd_target) )
  }
  
  conform_simulations[is.na(conform_simulations)] <- FALSE
  
  return(conform_simulations)
}

####################################

combined_likelihood <- function(conform_simulations, condition = 'all', marginal = 'neither'){
  
  if (condition == 'time_size_evenness'){
    combi_raw <- conform_simulations %>%
      rowwise() %>%
      mutate(all_conditions = time_CFU_recip & time_CFU_trans 
             & recip_pop_size & trans_pop_size
             & Even_recip & Even_trans)
  } else {
    return(NULL)
  }
  
  if (marginal == 'neither'){
    combi <- combi_raw %>%
      select(replicate, all_conditions, conj_rate, mig_rate) %>%
      group_by(conj_rate, mig_rate) %>%
      summarise(accepted_sims = sum(all_conditions),
                .groups = 'drop')
  } else if (marginal == 'raw'){
    combi <- combi_raw %>%
      select(replicate, all_conditions, conj_rate, mig_rate)
  } else if (marginal == 'mig_rate'){
    combi <- combi_raw %>%
      select(replicate, all_conditions, conj_rate, mig_rate) %>%
      group_by(mig_rate) %>%
      summarise(accepted_sims = sum(all_conditions),
                .groups = 'drop')
  } else if (marginal == 'conj_rate'){
    combi <- combi_raw %>%
      select(replicate, all_conditions, conj_rate, mig_rate) %>%
      group_by(conj_rate) %>%
      summarise(accepted_sims = sum(all_conditions),
                .groups = 'drop')
  }
  
  return(combi)
}
## Plotting functions #####
plot_combined_likelihood_function <- function(conform_simulations, sd_dist = 3, 
                                              condition = 'all'){
  
  combi <- combined_likelihood(conform_simulations, condition, marginal = 'neither')
  
  #Plot the result
  ggplot(data = combi, aes(x=conj_rate, y=mig_rate, fill= accepted_sims) )+ 
    geom_raster() +
    scale_x_continuous(trans='log10') + 
    scale_y_continuous(trans='log10') + 
    scale_fill_viridis(limits = c(0, 100)) + 
    annotation_logticks(base = 10, sides='lb') + 
    labs(x='Conjugation rate', 
         y='Migration rate', 
         fill='Percentage of \nsuccesful \nsimulations')
  
  # Save as png
  if (condition == 'all'){
    ggsave(paste0(plot_folder,"likeliness_sd_",sd_dist,".pdf"))
    write_csv(pivot_wider(combi, names_from = conj_rate, values_from = accepted_sims),
              paste0(simulation_data_folder, "likelihood_sd_",sd_dist,".csv"))
  } else {
    ggsave(paste0(plot_folder,"likeliness_sd_",sd_dist,"_", condition, ".pdf"))
    write_csv(pivot_wider(combi, names_from = conj_rate, values_from = accepted_sims),
              paste0(simulation_data_folder, "likelihood_sd_",
              sd_dist,"_", condition,".csv"))
  }
  
}

plot_all_likelihood_functions <- function(conform_simulations, sd_dist = 3){
  
  summary_stats = setdiff(colnames(conform_simulations), 
                          c("replicate", "conj_rate","mig_rate" ))
  
  for (stat in summary_stats){
  
  combi <- conform_simulations %>%
    select(replicate, condition = all_of(stat), conj_rate, mig_rate) %>%
    group_by(conj_rate, mig_rate) %>%
    summarise(accepted_sims = sum(condition),
              .groups = 'drop')
  
  #Plot the result
  ggplot(data = combi, aes(x=conj_rate, y=mig_rate, fill= accepted_sims) )+ 
    geom_raster() +
    scale_x_continuous(trans='log10') + 
    scale_y_continuous(trans='log10') + 
    scale_fill_viridis(limits = c(0, 100)) + 
    annotation_logticks(base = 10, sides='lb') + 
    labs(x='Conjugation rate', 
         y='Migration rate', 
         fill='Percentage of \nsuccesful \nsimulations') 
  
  # Save as png
  ggsave(paste0(plot_folder,"likeliness_sd_",sd_dist,"_", stat, ".pdf"))
  write_csv(pivot_wider(combi, names_from = conj_rate, values_from = accepted_sims),
              paste0(simulation_data_folder, "likelihood_sd_",sd_dist,"_", stat,".csv"))

  }
}

########################## Parameter Estimation ##############################################
conform_simulations_sd3 <- acceptance_condition(all_simulations, exp_stats, sd_dist = 3)

plot_combined_likelihood_function(conform_simulations_sd3, sd_dist = 3, condition = 'time_size_evenness')
plot_all_likelihood_functions(conform_simulations_sd3, sd_dist = 3)

####################################
# Estimates for the parameter values
library("HDInterval")

combi<-combined_likelihood(conform_simulations_sd3, condition = 'time_size_evenness')

# Max pair
combi[which.max(combi$accepted_sims),]

# HPD interval conj_rate ####
marginal_conj_rate <- combined_likelihood(conform_simulations_sd3, 
                                          condition = 'time_size_evenness', 
                                          marginal = 'conj_rate')

ggplot(marginal_conj_rate) +
  geom_line(aes(x=conj_rate, y = accepted_sims)) +
  scale_x_continuous(trans = 'log10')

test <- unlist(sapply(1:nrow(marginal_conj_rate), 
       function(x){rep(marginal_conj_rate[[x,'conj_rate']], 
                      marginal_conj_rate[x, 'accepted_sims'])}) )
hdi(test)
mean(test)

# HPD interval mig_rate ####

marginal_mig_rate <- combined_likelihood(conform_simulations_sd3, 
                                         condition = 'time_size_evenness',
                                         marginal = 'mig_rate')

ggplot(marginal_mig_rate) +
  geom_line(aes(x=mig_rate, y = accepted_sims)) +
  scale_x_continuous(trans = 'log10')

test2 <- unlist(sapply(1:nrow(marginal_mig_rate), 
                      function(x){rep(marginal_mig_rate[[x,'mig_rate']], 
                                      marginal_mig_rate[x, 'accepted_sims'])}) )
hdi(test2)
mean(test2)

# HPD interval both ####

marginal_raw <- combined_likelihood(conform_simulations_sd3, 
                                    condition = 'time_size_evenness',
                                    marginal = 'raw')
print(hdi(marginal_raw[marginal_raw$all_conditions == TRUE,]))


##################################



