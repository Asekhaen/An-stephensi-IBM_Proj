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

set.seed(20250606)

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
                   init_frequency = init_frequency,
                   bloodmeal_prob = bloodmeal_prob, 
                   fecundity = fecundity, 
                   conversion_prob,
                   resistance_prob,
                   daily_survival = daily_survival, 
                   daily_transition = daily_transition,
                   alpha = alpha,
                   beta = beta,
                   decay = decay,
                   fecundity_effect = fecundity_effect,
                   lethal_effect = FALSE,
                   complete_sterile = FALSE,
                   sim_days = sim_days,
                   dispersal_matrix = dispersal_matrix,
                   t_max,
                   t_min,
                   sigma,
                   gdd_required = gdd_required,
                   ldt = ldt)
 
