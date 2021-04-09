######################################################################################
## make_sim_plottable.R
##
## Transform the output of an adaptive tau simulation such that
## it can be plotted with ggplot (incl. subsetting for plasmid populations only)
##
## Author: Jana S. Huisman
## Last update: Feb 2021
######################################################################################

# Per plasmid tag population
make_tags_plottable = function(Ni, Nj, simulation){
  # Return the names of the plasmid populations in the simulated experiment
  populations = population_list(Ni, Nj)
  subset_plasmidtag <- pop_per_plasmid_tag(Nj, populations)
  
  # For each timestep, sum all populations with the same plasmid tag together
  if(length(subset_plasmidtag[[1]])>1){
    plasmid_pops <- lapply(1:Nj, function(i, simulation, subset_plasmidtag){rowSums(simulation[,subset_plasmidtag[[i]]])}, simulation, subset_plasmidtag)  
  }  else {
    plasmid_pops <- lapply(1:Nj, function(i, simulation, subset_plasmidtag){simulation[,subset_plasmidtag[[i]]]}, simulation, subset_plasmidtag)  
  }
  
  # Transform plasmid_pops into a data_frame with time column
  df <- as.data.frame(x= plasmid_pops, nrow=length(plasmid_pops[[1]]))
  colnames(df) <- rep(paste0("Plasmid ", 1:Nj))
  df<- cbind(simulation[,"time"], df)
  colnames(df)[1] <- "time"
  
  #Transform df into something plottable by ggplot
  df <- melt(df, id.vars="time")
  
  return(df)
}

