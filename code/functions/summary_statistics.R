######################################################################################
## summary_statistics.R
##
## Defines a set of measures to judge how close the simulations
## are to the experimental data
## Includes: Evenness, p0, rho, timing of the CFU trace
##
## Author: Jana S. Huisman
## Last update: Feb 2021
######################################################################################

CFU_timing <- function(simulation, population_names, marker_pop_size = 1e6){
  #sums all populations in population_names together
  summed_pops <- rowSums(simulation[,population_names])
  
  sim_pop_time <- simulation[which(summed_pops>=marker_pop_size)[1],"time"]
  return(sim_pop_time)
}

# Calculates the evenness score
evenness <- function(proportions, detection_limit = 8.9e-5){
  proportions[(proportions) < detection_limit] <- 0
  n <- length(proportions)
  if (sum(proportions)<=0){
    return(NA)
  } else{
    p <- sort(proportions/sum(proportions))
    
    totally.even <- sum(cumsum(rep(1/n, n))) #
    totally.uneven <- sum(c(rep(0, n-1),1))  #
    
    E <- (sum(cumsum(p)) - totally.uneven)/(totally.even - totally.uneven)
    return(E)
  }
}

# ratio between most and 2nd most abundant tag frequencies
tag_ratio <- function(proportions, detection_limit = 8.9e-5){
  proportions[(proportions)<detection_limit] <- 0
  sorted_p <- sort(proportions, decreasing = T)
  if (sorted_p[2] ==0){
    r <- NA
  } else {
    r <- as.double(sorted_p[1]/sorted_p[2])  
  }
  return(r)
}

###############################
# sum all populations per tag
pop_per_plasmid_tag <- function(simulation_row, Nj){
  populations <- names(simulation_row)
  
  subset_plasmidtag <- list()
  pop_plasmidtag <- vector(mode = 'numeric', length = Nj)
  for (j in 1:Nj){
      subset_plasmidtag[j] <- list(grep(pattern=paste0("T.", j, "$"), populations, value = T))
      pop_plasmidtag[j] <- sum(simulation_row[unlist(subset_plasmidtag[j])])
  }
  
  return(pop_plasmidtag)
}

pop_per_chromo_tag <- function(simulation_row, Ni){
  populations <- names(simulation_row)
  
  subset_chromotag <- list()
  pop_chromotag <- vector(mode = 'numeric', length = Ni)
  for (i in 1:Ni){
    subset_chromotag[i] <- list(c( grep(pattern=paste0("R", i), populations, value = T),
                                   grep(pattern=paste0("T", i, ".$"), populations, value = T) ))
    pop_chromotag[i] <- sum(simulation_row[unlist(subset_chromotag[i])])
  }
  
  return(pop_chromotag)
}

###############################
calculate_summary_statistics <- function(simulation, detection_limit = 8.9e-5){
  trans_population_names <- grep(pattern=paste0("T.."), colnames(simulation), value = T)
  recip_population_names <- grep(pattern=paste0("R."), colnames(simulation), value = T)
  donor_population_names <- grep(pattern=paste0("D."), colnames(simulation), value = T)
  
  Ni = length(recip_population_names)
  Nj = length(donor_population_names)
  
  CFU_trans <- CFU_timing(simulation, trans_population_names)
  CFU_recip <- CFU_timing(simulation, recip_population_names)
  
  #population sizes at the end of the simulation
  endpoint_ensemble_result <- simulation[dim(simulation)[1], 2:dim(simulation)[2]]
  trans_pop_size <- sum(endpoint_ensemble_result[trans_population_names])
  recip_pop_size <- sum(endpoint_ensemble_result[recip_population_names])
  
  ## Transconjugant quantities
  if (trans_pop_size==0){
    Even_trans <- NA
    prop_trans <- 0
    ratio_trans <- NA
    
  } else {
    # the relative proportion of plasmid tags
    plasmid_pops <- pop_per_plasmid_tag(endpoint_ensemble_result, Nj)
    t_proportions <- plasmid_pops/sum(plasmid_pops)
    t_proportions[(t_proportions) < detection_limit] <- 0
    
    Even_trans <- evenness(t_proportions)
    prop_trans <- sum((t_proportions) > detection_limit)/length(plasmid_pops)  
    ratio_trans <- tag_ratio(t_proportions)
  }
  
  ## Recipient quantities
  if (recip_pop_size==0 & trans_pop_size==0){
    Even_recip <- NA
    prop_recip <- 0
    ratio_recip <- NA
    
  } else {
    # the relative proportion of plasmid tags
    chromo_pops <- pop_per_chromo_tag(endpoint_ensemble_result, Ni)
    chromo_proportions <- chromo_pops/sum(chromo_pops)
    chromo_proportions[(chromo_proportions) < detection_limit] <- 0
    
    Even_recip <- evenness(chromo_proportions)
    prop_recip <- sum((chromo_proportions) > detection_limit)/length(chromo_pops)  
    ratio_recip <- tag_ratio(chromo_proportions)
  }
  
  return(c(CFU_trans, CFU_recip, Even_trans, prop_trans, ratio_trans, 
           Even_recip, prop_recip, ratio_recip, trans_pop_size, recip_pop_size))
}

###############################
## for use in the parallel ensemble simulation function
calculate_sim_measures <- function(i){
  detection_limit = 8.9e-5

  simulation = ssa.adaptivetau(initial_population, transitions, lvrates, params,
                               tf=time_max,
                               tl.params = list(epsilon=0.001))

  result <- calculate_summary_statistics(simulation, detection_limit)
  rm(simulation)

  return(result)
}

