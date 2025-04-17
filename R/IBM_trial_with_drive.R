library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse)

# TO DO!
# introduce a functional exponential dispersal model. DONE!
# implement density-dependence survival at the larval stage. DONE!
# introduce gene drive. ONGOING!
# introduce patch carrying capacity (larvae)
# introduce effect of temperature and humidity on transition/growth and survival

###########################################
#               PARAMETERS                #
###########################################

#set.seed(04042025)
set.seed(19042025)

# Call functions and parameters 
source("R/parameters.R")
source("R/IBM-functions.R")


###########################################
#            RUN   SIMULATION             #
###########################################

sim <- simulation (patches = patches,
                   n_per_patch = n_per_patch, 
                   coords = coords,
                   n_loci = n_loci,
                   mate_prob = mate_prob, 
                   bloodmeal_prob = bloodmeal_prob, 
                   fecundity = fecundity, 
                   drive_conversion_prob = drive_conversion_prob,
                   daily_survival = daily_survival, 
                   daily_transition = daily_transition,
                   carry_k = carry_k,
                   sim_days = sim_days,
                   dispersal_matrix = dispersal_matrix)
 
