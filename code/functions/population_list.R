######################################################################################
## population_list.R
##
## Creates a vector of population names needed for the adaptive tau simulations,
## as well as a vector of starting population sizes.
## Also contains some functions useful to subset results into chromosomal and plasmid populations
##
## Author: Jana S. Huisman
## Last update: Jan 2021
######################################################################################

# returns a list of population names, e.g. "D1" "R1" "T11"
population_list <- function(Ni, Nj){

  total_number_populations = Nj + Ni + Ni*Nj
  populations <- vector("character", total_number_populations)
  
  populations[1:Nj]<-paste0('D', 1:Nj)
  populations[(Nj+1):(Nj+Ni)]<- paste0('R', 1:Ni)
  
  counter <- Nj+Ni+1

  for (i in 1:Ni) {
    for (j in 1:Nj) {
      populations[[counter]]<-paste0('T', i, j)
      counter <- counter+1
    }
  }
  
  return(populations)
}

# Only recipient populations, i.e. Ri
recip_pop_names <- function(populations){
  return(grep(pattern=paste0("R.", "+$"), populations, value = T))
}

# Only transconjugant populations, i.e. Tij
trans_pop_names <- function(populations){
  return(grep(pattern=paste0("T.", "+$"), populations, value = T))
} 

# lists all populations per chromosomal tag
pop_per_chromo_tag <- function(Ni, populations){
  subset_chromosomaltag <- list()
  for (i in 1:Ni){
    subset_chromosomaltag[i] <- list(c(paste0("R", i),
                                     grep(pattern=paste0("T", i, ".$"), populations, value = T)))
  }
  return(subset_chromosomaltag)
}

# lists all populations per plasmid tag
pop_per_plasmid_tag <- function(Nj, populations, trans_only = TRUE){
  subset_plasmidtag <- list()
  for (j in 1:Nj){
    if (trans_only){
      subset_plasmidtag[j] <- list(grep(pattern=paste0("T.", j, "$"), populations, value = T))
    } else{
      subset_plasmidtag[j] <- list(c(paste0("D", j),
                                     grep(pattern=paste0("T.", j, "$"), populations, value = T)))
    }
    
  }
  return(subset_plasmidtag)
}


# To initialise the starting population size for each plasmid tag
generate_pWITS_pop_vector <-function(Nj, total_pop_size, split){
  
   pWITS_pop_vector <- vector("numeric", Nj)
  
  if (split =='equal'){
    pWITS_pop_vector[] = total_pop_size/Nj
  } else if(split =='1pct'){
    
    pWITS_pop_vector[] = total_pop_size/100.0
    pWITS_pop_vector[1] = (total_pop_size/100.0)*(100.0  - length(Nj)+1)
  }
  
   return(pWITS_pop_vector)  
}

# initialise all starting population sizes
init_population_list <- function(Ni, Nj, total_pop_size, split = 'equal'){
  
  population_names <- population_list(Ni, Nj)
  total_number_populations = Ni + Nj + Ni*Nj
  
  population_sizes <- vector(mode="numeric", total_number_populations)
  names(population_sizes) <- population_names
  
  population_sizes[1:Nj] <- generate_pWITS_pop_vector(Nj, total_pop_size, split)

  return(population_sizes)
}

