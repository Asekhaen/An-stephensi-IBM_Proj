# parameters

patches <- 10                        # Number of patches
n_per_patch <- c(1000, 1000, 1000, 1000, 1000, 0, 0, 0, 0, 0)                 # Initial number of individuals per patch
daily_survival <- c(egg = 0.8,      # daily survival prob
                    larva = 0.74, 
                    pupa = 0.73, 
                    adult = 0.80)    
daily_transition <- c(egg = 0.5,
                      larva = 0.5,
                      pupa = 0.5)  # daily transition prob
fecundity <- 10                      # Number of offspring per day per female mosquitoe 
carry_k <- 10000                     # Carrying capacity
mate_prob <- 0.75                   # Probability of mating
sim_days <- 100                     # Number of simulation in days
bloodmeal_prob <- 0.75              # Probability that a female find a blood meal

# dispersal parameters
lambda <- 5
dispersal_frac <- 0.05

# Gene Drive Parameters
n_loci <- 5
drive_conversion_prob <- 0.95  # Rate at which the drive allele converts the wild-type allele
