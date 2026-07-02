# install and load the following packages if needed

###########################################
#               PARAMETERS                #
###########################################

set.seed(246)

# Source functions and parameters 
source("R/dependencies.R")

# pop <- ini_pop(patches, n_per_patch, coords, n_loci, release_freq)


# plot(coords, cex = 4)
# text(coords, labels = 1:patches)


###########################################
#            RUN   SIMULATION             #
###########################################

output <- run_model (patches = patches,
                     n_per_patch = n_per_patch,
                     coords = coords,
                     n_loci = n_loci,
                     # stages = stages,
                     prob_wildtype1 = prob_wildtype1, 
                     prob_wildtype2 = prob_wildtype2,
                     per_release = per_release,
                     release_day = release_day,
                     bloodmeal_prob = bloodmeal_prob,
                     beta = beta,
                     decay = decay,
                     recomb = FALSE, # if TRUE = freely recombining, FALSE = no recombination
                     # lethal_effect = FALSE,
                     drive_effect = "sterility", # "sterility" or "x_shredding". Default is "no effect" 
                     drive_type = "homing",  #drive type: "homing", "toxin_antidote", or "x_shredder". Default is "Mendelian"
                     prob1 = cut_rate,
                     prob2 = homing_rate,
                     homing_prob = homing_prob,
                     shred_prob = shred_prob,
                     sim_days = sim_days,
                     dispersal_type = adjacency_matrix,  #neg_exponet_model or adjacency matix
                     t_max = temp_max,
                     t_min = temp_min,
                     sigma,
                     # surface_area = s_area,
                     dd_effect = dd_effect,
                     max_survival = max_survival,
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
#                                       beta = beta,
#                                       decay = decay,
#                                       fecundity_effect = fecundity_effect,
#                                       lethal_effect = FALSE,
#                                       complete_sterile = FALSE,
#                                       sim_days = sim_days,
#                                       dispersal_type = dispersal_matrix,
#                                       t_max,
#                                       t_min,
#                                       sigma,
#                                       gdd_required = gdd_required,
#                                       ldt = ldt)
#                                   }))






