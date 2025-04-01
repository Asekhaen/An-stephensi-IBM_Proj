# parameters

patches <- 6                        # Number of patches
n_per_patch <- c(500, 250, 50, 0, 0, 0)                 # Initial number of individuals per patch

daily_survival <- c(egg = 0.8,      # daily survival prob
                    larva = 0.54, 
                    pupa = 0.63, 
                    adult = 0.60)    
daily_transition <- c(egg = 0.5,
                      larva = 0.5,
                      pupa = 0.5)  # daily transition prob
fecundity <- 9                      # Number of offspring per day per female mosquitoe 
carry_k <- 1000                     # Carrying capacity
mate_prob <- 0.81                   # Probability of mating
sim_days <- 1000                      # Number of simulation in days
bloodmeal_prob <- 0.70              # Probability that a female find a blood meal
dispersal_frac <- 0.2