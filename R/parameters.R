# parameters

patches <- 5                              # Number of patches
n_per_patch <- c(1000,0,0,0,0)            # Initial number of individuals per patch
daily_survival <- c(egg = 0.8,            # daily survival prob
                    larva = 0.80, 
                    pupa = 0.75, 
                    adult = 0.80)  

daily_transition <- c(egg_larva = 0.8,
                      larva_pupa = 0.8,
                      pupa_adult = 0.8)   # daily transition prob

gdd_required <- c(egg_larva = 20,         # growth degree-days required for transition 
                  larva_pupa = 40,
                  pupa_adult = 30)

fecundity <- 50                           # Number of offspring per day per female mosquito
beta <- 50                                # the adult male population size at which the daily probability of mating is 0.5.
sim_days <- 150                            # Number of simulation in days
bloodmeal_prob <- 0.75                    # Probability that a female find a blood meal



# dispersal parameters
lambda <- 0.2
dispersal_frac <- 0.05



# Genetic (load) & drive parameters
n_loci <- 10
init_frequency = 0.1                    # initial frequency of deleterious recessives
conversion_prob <- 0.95                  # Rate at which the drive allele converts the wild-type allele
resistance_prob <- 0.5                   # prob resistance development or conversion failure
fecundity_effect <- 1                    # effect per homozygous deleterious recessive on fecundity. Set to '0' to turn off
decay <- 0.5                             


# Environmental/ecological


ldt <- 14

t_max  <-  matrix(rnorm(patches * sim_days, 
                       mean = 30, sd = 5), 
                       nrow = sim_days, ncol = patches)  

t_min <-  matrix(rnorm(patches * sim_days, 
                       mean = 20, sd = 4), 
                 nrow = sim_days, ncol = patches)

  
sigma <- 8                              # Decay parameter for temperature effect on survival

alpha <- 0.0001                         # strength of density dependence 
