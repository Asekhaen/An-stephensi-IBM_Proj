# parameters

patches <- 3                         # Number of patches
n_per_patch <- c(5000, 0, 0)         # Initial number of individuals per patch
daily_survival <- c(egg = 0.8,       # daily survival prob
                    larva = 0.74, 
                    pupa = 0.73, 
                    adult = 0.80)    
daily_transition <- c(egg = 0.5,
                      larva = 0.5,
                      pupa = 0.5)   # daily transition prob
fecundity <- 10                     # Number of offspring per day per female mosquitoe 
beta <- 100                         # assumed to be the male population at which the probability of mating = 0.5 (North and Godfray Malar J (2018) 17:140)
sim_days <- 100                     # Number of simulation in days
bloodmeal_prob <- 0.70              # Probability that a female find a blood meal

# dispersal parameters
lambda <- 0.2
dispersal_frac <- 0.05

# Gene Drive Parameters
n_loci <- 5
conversion_prob <- 0.95  # Rate at which the drive allele converts the wild-type allele
resistance_prob <- 0.5   # prob resistance development or conversion failure


# Environmental/ecological
carry_k <- 50000                     # Carrying capacity
temp  <-  round((rnorm(sim_days, 
                       mean = 25, 
                       sd = 3)), 1)
sigma <- 8


