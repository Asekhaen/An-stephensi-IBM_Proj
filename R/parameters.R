

# general parameters ----------------------------------------------------------

set.seed(123)

# simulation time steps
sim_days <- 20 

#################################################
# life history parameters
#################################################

# Number of patches
patches <- 10     

# Initial number of individuals per patch
initial_pop = 3000                  # carrying capacity
carrying_capacity = 10000                  # carrying capacity
n_per_patch <- create_n_per_patch(patches, initial_pop)   # Initial number of individuals per patch
sex_alleles <- c("X", "Y") 
stages <- c("egg", "larva", "pupa", "adult")


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

dispersal_prop <- 0.01


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



cut_rate = 0.75      #cleavage rates: 50 ~ 90% (75) (the others 25% follow  Mendelian inheritance pattern)  
conv_eff = runif(1, min = 50, max = 95)    #conversion efficiency rates = 50 ~ 95% (75)
r1_eff = 0.35         #In-frame resistance development rate (restore Mendelian inheritance) = 10% 
r2_eff = 0.15         #Out-frame resistance dev (fitness is lost and drive can’t cut) = 6%




n_loci <- 2   # number of loci
init_frequency = 0.25     # initial frequency of deleterious recessives
cleavage = 0.95     # gene cleavage probability 
homing_rate =  0.95
conversion_prob <- cleavage * homing_rate    # Rate at which the drive allele converts the wild-type allele
resistance_prob <- 0.5   # prob resistance development or conversion failure
fecundity_effect <- 0   # effect per homozygous deleterious recessive on fecundity. 0 = no effect on fecundity or batch size. > 0 = 1 additive effect. 
decay <- 0.5     #controls the rate at which the covariance between two loci decreases with distance
l.cov.mat <- place_loci_mat(n_loci, genome.size = 1, var = 1, decay)


################################################
# environmental parameters
################################################

# create coordinates and dipersal matrix for the patches/locations 

# random location coordinates
set.seed(234)
coords <- as.data.frame(100 * matrix(runif(patches * 2), ncol = 2))
colnames(coords) <- c("x","y")
plot(coords)

# dispersal matrix
neg_exponet_model <- metapop (coords = coords, 
                                 lambda = lambda, 
                                 disp_prob = dispersal_prop)

adjacency_matrix <- step_stone(n_patches = patches, 
                               disp_prob = dispersal_prop)

# This bit of code generates random daily temperature and humidity to estimates  
# growth degree-day required for each stage during transition and survival
temp_max  <-  matrix(runif(patches * sim_days, 
                        min = 28, max = 38), 
                  nrow = sim_days, ncol = patches) 

temp_min <-  matrix(runif(patches * sim_days, 
                       min = 19, max = 25), 
                 nrow = sim_days, ncol = patches)

humidity <- matrix(rtruncnorm(patches * sim_days, a = 0, b = 100,
                        mean = 80, sd = 10), 
                  nrow = sim_days, ncol = patches)

# a coefficient that controls the strength of density dependence 
alpha <- 0.0001   

# size of the habitat in cm*
s_area <- 10000                   





