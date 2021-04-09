###########################################
## lvrates.R
##
## Returns a function that encodes the transition
## rates of migration, birth, death, and plasmid transmission
## processes for use in the adaptive-tau simulation
##
## This function includes a carrying capacity.
##
## Author: Jana S. Huisman
## Last update: Jan 2021
###########################################

# When specifying an experiment with 2 chromosomal tags
# and 1 tagged plasmid; i.e., Ni=2, Nj=1
# lvrates_carrying_cap_function(2,1) 
# the resulting function looks as follows:

# function(x, params, t) {
#   return(c(x["D1"]*params$birth_rate*(1+0.0833333333333333), 
#            x["R1"]*params$birth_rate*(1+0.0833333333333333) + params$migration_rate, 
#            x["R2"]*params$birth_rate*(1+0.0833333333333333) + params$migration_rate, 
#            x["T11"]*params$birth_rate*(1+0.0833333333333333), 
#            x["T21"]*params$birth_rate*(1+0.0833333333333333), 
#            x["D1"]*(params$death_rate+((x["D1"])/1e+09)*params$birth_rate), 
#            x["R1"]*(params$death_rate+((x["R1"]+x["R2"]+x["T11"]+x["T21"])/1e+09)*params$birth_rate), 
#            x["R2"]*(params$death_rate+((x["R1"]+x["R2"]+x["T11"]+x["T21"])/1e+09)*params$birth_rate), 
#            x["T11"]*(params$death_rate+((x["R1"]+x["R2"]+x["T11"]+x["T21"])/1e+09)*params$birth_rate), 
#            x["T21"]*(params$death_rate+((x["R1"]+x["R2"]+x["T11"]+x["T21"])/1e+09)*params$birth_rate), 
#            (params$conj_donor*x["D1"]+params$conj_trans*(x["T11"]+x["T21"]))*x["R1"],
#            (params$conj_donor*x["D1"]+params$conj_trans*(x["T11"]+x["T21"]))*x["R2"]))     }

#The order of these rates must match the order 
# in which transitions were specified.
#D, R, T - birth, death, transition

###########################################

lvrates_function <- function(Ni, Nj, carrying_cap_D = 1e9, 
                             carrying_cap_R = 1e9, leftover_birth = 1./12){
  chromosomal_range = 1:Ni
  plasmid_range = 1:Nj
  
  ######################
  # To pre-write some lengthy terms needed for the 
  # carrying capacity term in the Death rates
  all_d_population_terms <- paste0('(', paste0('x["D',1:Nj,'"]', collapse = '+'), ')' )
  #all_r_population_terms <- paste0('(', paste0('x["R',1:Ni,'"]', collapse = '+'), ')' )
  
  # so that there is an initial content that we add to (in the loop when i=1)
  all_rt_population_terms <- paste0('(', paste0('x["R',1:Ni,'"]', collapse = '+') )
  
  # j = 1, 2, ..., Nj
  for (j in plasmid_range){
    for (i in chromosomal_range) {
      all_rt_population_terms <- paste0(all_rt_population_terms, '+x["T',i,j,'"]')
    }
  }
  all_rt_population_terms <- paste0(all_rt_population_terms, ')')
  
  ######################
  function_def <- 'function(x, params, t) {
  return(c('
  
  ######################
  # Birth processes
  function_def <- paste0(function_def, 
                        paste0('x["D', 1:Nj, '"]*params$birth_rate*(1+',leftover_birth,'), ', collapse = ''))
  
  function_def <- paste0(function_def, 
                         paste0('x["R', 1:Ni, '"]*params$birth_rate*(1+',leftover_birth,') + params$migration_rate/', Ni, ', ', collapse = ''))
  
  for (i in chromosomal_range) {
    for (j in plasmid_range) {
        function_def<-paste0(function_def, 'x["T', i, j, '"]*params$birth_rate*(1+',leftover_birth,'), ')  
    }
  }
  ######################
  # Death processes
  function_def <- paste0(function_def, 
                         paste0('x["D', 1:Nj, '"]*(params$death_rate+(',all_d_population_terms,'/',carrying_cap_D,')*params$birth_rate), ', collapse = ''))
  
  function_def <- paste0(function_def, 
                         paste0('x["R', 1:Ni, '"]*(params$death_rate+(',all_rt_population_terms,'/',carrying_cap_R,')*params$birth_rate), ', collapse = ''))
  
  for (i in chromosomal_range) {
    for (j in plasmid_range) {
      function_def<-paste0(function_def, 'x["T', i, j, '"]*(params$death_rate+(',all_rt_population_terms,'/',carrying_cap_R,')*params$birth_rate), ')
    }
  }
  ######################
  # To pre-write some lengthy terms
  # each list item contains a string with all transconjugant populations
  # with that plasmid
  transconjugant_terms <- lapply(1:Nj, function(j) {paste0('x["T',1:Ni,j,'"]', collapse = '+')})
  ######################
  #Interaction
  for (i in chromosomal_range) {
    for (j in plasmid_range) {
      function_def<-paste0(function_def, '(params$conj_donor*x["D',j,'"]',
                           '+params$conj_trans*(', transconjugant_terms[j], '))*x["R',i,'"]', ',')
    }
  }
  # Because the last entry should not end with a ',' but with function end instead
  function_end <- '))     }'
  function_def <- sub(pattern=',$', replacement=function_end, x=function_def)
  
  return(eval(parse(text=function_def)))
}
