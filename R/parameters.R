# parameters

patches <- 6                              # Number of patches
n_per_patch <- c(5000,0,0,0,0,0)            # Initial number of individuals per patch
daily_survival <- c(egg = 0.8,            # daily survival prob
                    larva = 0.80,
                    pupa = 0.80,
                    adult = 0.80)

# daily_transition <- c(egg_larva = 0.8,
#                       larva_pupa = 0.8,
#                       pupa_adult = 0.8)   # daily transition prob


fecundity <- 10                           # Number of offspring per day per female mosquito
beta <- 100                                # the adult male population size at which the daily probability of mating is 0.5.
sim_days <- 100                           # Number of simulation in days
bloodmeal_prob <- 0.85                    # Probability that a female find a blood meal



# dispersal parameters
lambda <- 0.1
dispersal_frac <- 0.01



# Genetic (load) & drive parameters
n_loci <- 5
init_frequency = 0.2                    # initial frequency of deleterious recessives
conversion_prob <- 0.95                  # Rate at which the drive allele converts the wild-type allele
resistance_prob <- 0.5                   # prob resistance development or conversion failure
fecundity_effect <- 0.2             # effect per homozygous deleterious recessive on fecundity. Set to '0' to turn off
decay <- 0.5                             


# Environmental/ecological


ldt <- c(egg = 8.19, 
         larva = 9.33, 
         pupa = 10.3,
         adult = 11.5)

t_max  <-  matrix(rnorm(patches * sim_days, 
                       mean = 30, sd = 5), 
                       nrow = sim_days, ncol = patches)  

t_min <-  matrix(rnorm(patches * sim_days, 
                       mean = 20, sd = 5), 
                 nrow = sim_days, ncol = patches)
sigma <- 8                              # Decay parameter for temperature effect on survival
alpha <- 0.0001                         # a coefficient that controls the strength of density dependence 
surface_area <- 10000                   # size of the habitat*
