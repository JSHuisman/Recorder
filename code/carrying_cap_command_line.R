###########################################################
## RECORDER - sweep across the carrying capacity
##
## Author: J.S. Huisman
## Last update: Feb 2021
###########################################################
# To get results for higher/lower conj rate, change the
# values in line 38 to 42/43, and save elsewhere (line 76)
###########################################################
library("adaptivetau")

source("functions/transition_list.R")
source("functions/lvrates.R")
source("functions/population_list.R")

source("functions/summary_statistics.R")
source("functions/parallel_ensemble_sim_function.R")

#################################################################
input_args <- commandArgs()
# exclude R directory and R commands
condition = !(grepl(pattern="^--", input_args)|grepl(pattern="^/", input_args))
user_inputs <- input_args[condition]

number_of_cores = as.numeric(user_inputs[1])
carrying_cap_index = as.numeric(user_inputs[2]) 

experiment=paste0(user_inputs[2])
#################################################################
#Parameters for all simulations

birth_rate = 44*log(2)
death_rate = 4*log(2)
leftover_birth = 1./12
total_init_pop = 1e7

# Best estimates
conj_rate = 10**(-10.5)
mig_rate = 1.78

# Additional values [higher, lower]
#conj_rate = 10**(-8.5)
#conj_rate = 10**(-12.5)

params = list(conj_donor= conj_rate, 
              conj_trans= conj_rate, 
              birth_rate= birth_rate, death_rate= death_rate,
              migration_rate= mig_rate)

Ni = 3 #chromosomal tags
Nj = 4 #plasmid tags

#################################################################
ensemble_size = 100
carrying_cap_options = 10**seq(3,10,0.25)
#################################################################
#Setting up the simulations

transitions = transition_list(Ni, Nj)
lvrates = lvrates_function(Ni, Nj, 
                           carrying_cap_D = carrying_cap_options[carrying_cap_index], 
                           carrying_cap_R = carrying_cap_options[carrying_cap_index], 
                           leftover_birth)
populations = population_list(Ni, Nj)

# starting values
pWITS_distribution_type = 'equal' 
initial_population = init_population_list(Ni, Nj, total_init_pop, pWITS_distribution_type)

time_steps = 1000
time_max = 8
##################################################################
result<- parallel_ensemble_sim_function(ensemble_size, initial_population, transitions, lvrates,
                                        params,number_of_cores)
write.csv(result, paste0('../simulations_cc/', experiment, '.csv'))
#write.csv(result, paste0('../simulations_cc_high_conj/', experiment, '.csv'))
#write.csv(result, paste0('../simulations_cc_low_conj/', experiment, '.csv'))