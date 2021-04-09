######################################################################################
## parallel_ensemble_sim_function.R
##
## Calculates the summary statistics for an ensemble of simulations
##
## Author: Jana S. Huisman
## Last update: Feb 2021
######################################################################################

# requires summary_statistics.R to be loaded to find calculate_sim_measures()

parallel_ensemble_sim_function <- function(ensemble_size, initial_population, transitions, lvrates,
                                           params, number_of_cores = 24){
  library("parallel")

  # SUBSTITUTE this by your own number of cores, and/or local computer!
  if (as.logical(Sys.info()[1] == "Darwin")){
    num.cores <- max(1, detectCores() - 1) 
  } else {
    num.cores <- number_of_cores - 1
  }
  
  # Initialise the cluster
  cl <- makeCluster(num.cores, type="FORK")
  # Set a seed for the random number generator, which then seeds all new jobs
  # to ensure no random sequences are used on multiple cores
  clusterSetRNGStream(cl, 150)
  clusterExport(cl=cl, list(initial_population="initial_population", 
                            transitions="transitions", lvrates="lvrates",
                            params="params"))
  
  # Now our function
  result <- parSapply(cl, seq_len(ensemble_size), calculate_sim_measures, simplify = 'matrix')
  
  row.names(result)<- c("time_CFU_trans", "time_CFU_recip", "Even_trans", "prop_trans", "ratio_trans", 
                        "Even_recip", "prop_recip", "ratio_recip", "trans_pop_size", "recip_pop_size")

  # Finally, close the cluster we created
  stopCluster(cl)  
  # clean up memory because otherwise cluster connection problems can occur
  gc()
  return(t(result))
}
