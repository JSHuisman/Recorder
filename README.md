# Pathogen invasion-dependent tissue reservoirs and plasmid-encoded antibiotic degradation boost plasmid spread in the gut

This repository contains the data and code accompanying the manuscript "Pathogen invasion-dependent tissue reservoirs and plasmid-encoded antibiotic degradation boost plasmid spread in the gut" by Erik Bakkeren, Joana A. Herter, Jana S. Huisman, Yves Steiger, Ersin Gül, Joshua P.M. Newson, Alexander O. Brachmann, Jörn Piel, Roland R. Regoes, Sebastian Bonhoeffer, Médéric Diard, and Wolf-Dietrich Hardt.

## Experimental Data
The folder "source_data" contains excel sheets with all data needed to plot the graphs, and includes exact p values for statistical tests.
The folder "data" contains the data in a format compatible with the simulation/likelihood calculation code.

## Code
We have set up stochastic simulations to model the population dynamics of Salmonella in the murine gut, with reseeding from systemic sites. The model is primarily defined in the scripts in the functions folder (see explanation below). 

*Functions*
- lvrates_carrying_cap.R: Returns a function that encodes the transition rates of birth, death, and plasmid transmission processes for use in the adaptive-tau simulation.
- population_list.R: Creates a vector of populations needed for the adaptive tau simulations, as well as a vector of starting population sizes. Also contains some functions useful to subset results into chromosomal and plasmid populations
- transition_list.R: Creates a list of birth, death, and plasmid transmission transitions for use in the adaptive-tau simulation.
- summary_statistics.R: Defines a set of measures to judge how close the simulations are to the experimental data. 
- parallel_ensemble_sim_function.R: Calculates a set of measures from an ensemble of simulations
- make_sim_plottable.R: Transforms the output of a single adaptive tau simulation such that it can be plotted with ggplot.

**Note: Please set the maximum number of cores you wish to use in the parallel computation in functions/parallel_ensemble_sim_function.R**

## Simulation Results
Since the simulations take quite some computing power to reproduce, we have added the results in the folders called simulations*.





