# parameters

patches <- 5                         # Number of patches
n_per_patch <- c(5000,0,0,0,0)         # Initial number of individuals per patch
daily_survival <- c(egg = 0.8,       # daily survival prob
                    larva = 0.80, 
                    pupa = 0.75, 
                    adult = 0.80)    
daily_transition <- c(egg_larva = 0.5,
                      larva_pupa = 0.5,
                      pupa_adult = 0.5)    # daily transition prob
fecundity <- 10                      # Number of offspring per day per female mosquito
beta <- 50                           # the adult male population size at which the daily probability of mating is 0.5.
sim_days <- 50                       # Number of simulation in days
bloodmeal_prob <- 0.75               # Probability that a female find a blood meal

# dispersal parameters
lambda <- 0.2
dispersal_frac <- 0.05

# Gene Drive Parameters
n_loci <- 10
init_frequency = 0.75 # initial frequency of deleterious recessives
conversion_prob <- 0.95  # Rate at which the drive allele converts the wild-type allele
resistance_prob <- 0.5   # prob resistance development or conversion failure
fecundity_effect <- 0.2 # effect per homozygous deleterious recessive on fecundity. Set to zero to turn off
decay <- 0.5


# Environmental/ecological
temp  <-  round((rnorm(sim_days, 
                       mean = 25, 
                       sd = 0)), 1)
sigma <- 8        # Decay parameter for temperature effect on survival

alpha <- 0.0001   # strength of density dependence 
