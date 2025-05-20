library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse)

# TO DO!
# introduce a functional exponential dispersal model. DONE!
# implement density-dependence survival at the larval stage. DONE!
# introduce patch carrying capacity (larvae). DONE!
# introduce effect of temperature on survival (larva and adult) DONE!
# introduce effect of temperature and humidity on transition/growth
# introduce gene drive. ONGOING

###########################################
#               PARAMETERS                #
###########################################

set.seed(20250512)

# Source functions and parameters 
source("R/parameters.R")
source("R/IBM_functions.R")

plot(coords, cex = 4)
text(coords, labels = 1:patches)


###########################################
#            RUN   SIMULATION             #
###########################################

sim <- simulation (patches = patches,
                   n_per_patch = n_per_patch, 
                   coords = coords,
                   n_loci = n_loci, 
                   bloodmeal_prob = bloodmeal_prob, 
                   fecundity = fecundity, 
                   conversion_prob,
                   resistance_prob,
                   daily_survival = daily_survival, 
                   daily_transition = daily_transition,
                   beta = beta,
                   alpha = alpha,
                   sim_days = sim_days,
                   dispersal_matrix = dispersal_matrix,
                   daily_temp = temp[day],
                   sigma)
 
