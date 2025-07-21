library(tidyverse)
library(purrr)
library(truncnorm)


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

set.seed(20250711)

# Source functions and parameters 
source("R/parameters.R")
source("R/IBM_functions.R")
#source("R/IBM_functions_gdd.R")


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
                   gdd_required = gdd_required,
                   lethal_effect = TRUE,
                   complete_sterile = FALSE,
                   sim_days = sim_days,
                   stepping_stone_model = FALSE,
                   dispersal_matrix = dispersal_matrix,
                   t_max,
                   t_min,
                   sigma,
                   surface_area = surface_area,
                   ldt = ldt,
                   mu = mu,
                   sigma_dd = sigma_dd)
 


# # To run multiple scenarios with varying parameters, use the #purrr:pmap" function
# 
# 
# # Example: different init_frequency and fecundity
# simulation_scenarios <- expand.grid(
#   init_frequency = c(0.05, 0.1, 0.25, 0.5),
#   fecundity = c(5, 10, 20, 50),
#   fecundity_effect = c(0, 0.1, 0.2)
#   n_loci = c(2,4,8,16,32)
# )
# 
# 
# sim_output <- simulation_scenarios |>
#   mutate(simulation_result = pmap(.l = list(init_frequency, fecundity, fecundity_effect, n_loci),
#                                   .f = function(init_frequency, fecundity, fecundity_effect, n_loci) {
#                                     simulation (
#                                       patches = patches,
#                                       n_per_patch = n_per_patch, 
#                                       coords = coords,
#                                       n_loci = n_loci, 
#                                       init_frequency = init_frequency,
#                                       bloodmeal_prob = bloodmeal_prob, 
#                                       fecundity = fecundity, 
#                                       conversion_prob,
#                                       resistance_prob,
#                                       daily_survival = daily_survival, 
#                                       daily_transition = daily_transition,
#                                       alpha = alpha,
#                                       beta = beta,
#                                       decay = decay,
#                                       fecundity_effect = fecundity_effect,
#                                       lethal_effect = FALSE,
#                                       complete_sterile = FALSE,
#                                       sim_days = sim_days,
#                                       dispersal_matrix = dispersal_matrix,
#                                       t_max,
#                                       t_min,
#                                       sigma,
#                                       gdd_required = gdd_required,
#                                       ldt = ldt)
#                                   }))






