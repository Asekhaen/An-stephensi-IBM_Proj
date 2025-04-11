library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse)

# TO DO!
# introduce a functional exponential dispersal model. DONE!
# implement density-dependence survival at the larval stage. DONE!
# introduce gene drive. DONE!
# introduce patch carrying capacity 
# introduce effect of temperature and humidity on transition/growth and survival

###########################################
#               PARAMETERS                #
###########################################

#set.seed(04042025)
set.seed(04042022)

source("R/IBM-functions.R")
source("R/parameters.R")


###########################################
#  INITIALISE POPULATION WITH ATTRIBUTES  #
###########################################

# create coordinates for the patches/locations 
coords <- as.data.frame(100 * matrix(runif(patches * 2), ncol = 2))
colnames(coords) <- c("x","y")

# create a dispersal matrix
dispersal_matrix <- make_dispersal_matrix(coords = coords, 
                                          lambda = lambda, 
                                          dispersal_frac = dispersal_frac)

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
 
