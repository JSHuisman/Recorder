###########################################
## transition_list.R
##
## Creates a list of birth, death, and plasmid transmission
## transitions for use in the adaptive-tau simulation
##
## Author: Jana S. Huisman
## Last update: Jan 2021
###########################################
# the resulting list looks as follows:
#
# list(c(D1 = 1), # D1 birth
#      c(D2 = 1), # D2 birth
#      c(R1 = 1) # R1 birth /migration
#      c(T11 = 1), # T11 birth
#      c(T12 = 1), # T12 birth
#      c(D1 = -1), # D1 death
#      c(D2 = -1), # D2 death
#      c(R1 = -1) # R1 death
#      c(T11 = -1), # T11 death
#      c(T12 = -1), # T12 death
#      c(R1 = -1, T11 = +1), # plasmid donation from donor or other transconjugant
#      c(R1 = -1, T12 = +1), # plasmid donation from donor or other transconjugant
# )

transition_list <- function(Ni, Nj){
  #population names, e.g. "D1" "R1" "T11"
  
  i_range = 1:Ni
  j_range = 1:Nj #plasmid
  
  n_pop = Nj + Ni + Ni*Nj
  total_number_transitions = 2*n_pop + Ni*Nj #birth, death, transfer R->T
  
  transitions <- vector("list", total_number_transitions)
  
  transitions[1:Nj]<- lapply(1:Nj, function(x){setNames(1, paste0('D', x) )} )
  transitions[(Nj+1):(Nj+Ni)]<- lapply(1:Ni, function(x){setNames(1, paste0('R', x) )} )
  
  
  counter <- Nj+Ni +1
  # Birth
  for (i in i_range) {
    for (j in j_range) {
      transitions[[counter]]<-eval(parse(text=paste0('c(T', i, j, '= +1)')))
      counter <- counter+1
    }
  }
  
  #Death
  transitions[counter:(counter+Nj-1)]<- lapply(1:Nj, function(x){setNames(-1, paste0('D', x) )} )
  transitions[(counter+Nj):(counter+Nj+Ni-1)]<- lapply(1:Ni, function(x){setNames(-1, paste0('R', x) )} )
  
  counter <- counter + Nj+Ni
  for (i in i_range) {
    for (j in j_range) {
      transitions[[counter]]<-eval(parse(text=paste0('c(T', i, j, '= -1)')))
      counter <- counter+1
    }
  }
  
  #Interaction
  for (i in i_range) {
    for (j in j_range) {
      transitions[[counter]]<-eval(parse(text=paste0('c(R', i,'= -1, T', i, j, '= +1)')))
      counter <- counter+1
    }
  }
  
  return(transitions)
}

