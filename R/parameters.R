
# parameters

# simulation time steps
sim_days <- 25 

#################################################
# life history parameters
#################################################

# Number of patches
patches <- 10     

# Initial number of individuals per patch
n_per_patch <- c(5000,0,0,0,0,0,0,0,0,0)   

# the adult male population size at which the daily probability of mating is 
# 0.5 (North and Godfray, Malar J (2018) 17:140)
beta <- 100      

# Probability that a female find a blood meal
bloodmeal_prob <- 0.40                   

# low degree-day threshold
ldt <- c(egg = 8.19, 
         larva = 9.33, 
         pupa = 10.3)

# dispersal parameters

# lambda controls the rates at which probability between patches decreases with distances 
lambda <- 0.1

dispersal_prop <- 0.002



################################################
# Genetic (load) & drive parameters
#################################################


# number of loci
n_loci <- 5

# initial frequency of deleterious recessives
init_frequency = 0.25                   

# Rate at which the drive allele converts the wild-type allele
conversion_prob <- 0.95       

# prob resistance development or conversion failure
resistance_prob <- 0.5                   

# effect per homozygous deleterious recessive on fecundity. 0 = no effect on fecundity or batch size. > 0 = 1 additive effect. 
fecundity_effect <- 0         

#controls the rate at which the covariance between two loci decreases with distance
decay <- 0.5                            



################################################
# environmental parameters
################################################

# create coordinates for the patches/locations 
coords <- as.data.frame(100 * matrix(runif(patches * 2), ncol = 2))
colnames(coords) <- c("x","y")


# saved model used to estimate Life history parameters (Temperature and humidity) 
# for vector survival estimated from Golding et al., unpublished)

aquatic_stage <- "C:/Users/22181916/Documents/Curtin-PhD/R_and_IBM/An-stephensi-IBM_Proj/R/das_temp_dens_As.RDS"
adult_stage <- "C:/Users/22181916/Documents/Curtin-PhD/R_and_IBM/An-stephensi-IBM_Proj/R/ds_temp_humid.RDS"


# This bit of code generates random daily temperature and humidity to estimates  
# growth degree-day required for each stage during transition and survival
t_max  <-  matrix(rnorm(patches * sim_days, 
                       mean = 30, sd = 5), 
                       nrow = sim_days, ncol = patches)  

t_min <-  matrix(rnorm(patches * sim_days, 
                       mean = 20, sd = 5), 
                 nrow = sim_days, ncol = patches)

humidity <- matrix(rtruncnorm(patches * sim_days, a = 0, b = 100,
                        mean = 80, sd = 10), 
                  nrow = sim_days, ncol = patches)

# a coefficient that controls the strength of density dependence 
alpha <- 0.0001   

# size of the habitat in cm*
surface_area <- 10000                   


