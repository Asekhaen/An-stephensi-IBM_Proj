

# general parameters ----------------------------------------------------------



# simulation time steps
sim_days <- 50 

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


# growth degree day parameters
# stage development using degree-days (dd): degree-days, percentage/probability of 
# transitioning (data from Abbasi et al., 2023)   

# Eggs transition
one_per_egg <- sigma_etimate(30.5, 44.5, 0.01)
ten_per_egg <- sigma_etimate(33.7, 44.5, 0.1)
eighty_per_egg <- sigma_etimate(59.6, 44.5, 0.8)
mean_sigma_egg <- mean(c(one_per_egg, ten_per_egg, eighty_per_egg))

#Larva transition
one_per_larva <- sigma_etimate(93.9, 145.5, 0.01)
ten_per_larva <- sigma_etimate(104.1, 145.5, 0.1)
eighty_per_larva <- sigma_etimate(182.3, 145.5, 0.8)
mean_sigma_larva <- mean(c(one_per_larva, ten_per_larva, eighty_per_larva))

#Pupa transition
one_per_pupa <- sigma_etimate(17.7, 29.7, 0.01)
ten_per_pupa <- sigma_etimate(22.3, 29.7, 0.1)
eighty_per_pupa <- sigma_etimate(40.4, 29.7, 0.8)
mean_sigma_pupa <- mean(c(one_per_pupa, ten_per_pupa, eighty_per_pupa))


mu <- c(egg = 44.5, larva = 145.5, pupa = 29.7)
sigma_dd <- c(egg = mean_sigma_egg, 
              larva = mean_sigma_larva, 
              pupa = mean_sigma_pupa)


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

l.cov.mat <- place_loci_mat(n_loci, genome.size = 1, var = 1, decay)


################################################
# environmental parameters
################################################

# create coordinates and dipersal matrix for the patches/locations 

# random location coordinates
coords <- as.data.frame(100 * matrix(runif(patches * 2), ncol = 2))
colnames(coords) <- c("x","y")


# dispersal matrix
dispersal_matrix <- make_dispersal_matrix(coords = coords, 
                                          lambda = lambda, 
                                          dispersal_prop = dispersal_prop)


# This bit of code generates random daily temperature and humidity to estimates  
# growth degree-day required for each stage during transition and survival
t_max  <-  matrix(rnorm(patches * sim_days, 
                       mean = 30, sd = 5), 
                       nrow = sim_days, ncol = patches)  

t_min <-  matrix(rnorm(patches * sim_days, 
                       mean = 28, sd = 3), 
                 nrow = sim_days, ncol = patches)

humidity <- matrix(rtruncnorm(patches * sim_days, a = 0, b = 100,
                        mean = 80, sd = 10), 
                  nrow = sim_days, ncol = patches)

# a coefficient that controls the strength of density dependence 
alpha <- 0.0001   

# size of the habitat in cm*
s_area <- 10000                   





