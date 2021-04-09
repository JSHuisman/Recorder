###########################################################
## RECORDER - on EULER
##
## Author: J.S. Huisman
## Last update: Feb 2021
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
conj_rate_index = as.numeric(user_inputs[2]) 
mig_rate_index = as.numeric(user_inputs[3]) 

experiment=paste0(user_inputs[2],"_", user_inputs[3])
#################################################################
#Parameters for all simulations

birth_rate = 44*log(2)
death_rate = 4*log(2)
carrying_cap = 1e9
leftover_birth = 1./12
total_init_pop = 1e7

Ni = 3 #chromosomal tags
Nj = 4 #plasmid tags

#################################################################
#Setting up the simulations

transitions = transition_list(Ni, Nj)
lvrates = lvrates_function(Ni, Nj, carrying_cap_D = carrying_cap, 
                           carrying_cap_R = carrying_cap, leftover_birth)
populations = population_list(Ni, Nj)

# starting values
pWITS_distribution_type = 'equal' 
initial_population = init_population_list(Ni, Nj, total_init_pop, pWITS_distribution_type)

time_steps = 1000
time_max = 8
#################################################################
ensemble_size = 100
#the length of these vectors determines how to loop in the 
# accompanying shell script
conj_rate_options = 10**seq(-12,-3,0.25) #37
mig_rate_options = 10**seq(-6,3,0.25) #37

params = list(conj_donor= conj_rate_options[conj_rate_index], 
              conj_trans= conj_rate_options[conj_rate_index], 
              birth_rate= birth_rate, death_rate= death_rate,
              migration_rate= mig_rate_options[mig_rate_index])

##################################################################
result<- parallel_ensemble_sim_function(ensemble_size, initial_population, transitions, lvrates,
                                             params,number_of_cores)
write.csv(result, paste0('../simulations/', experiment, '.csv'))