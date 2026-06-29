

# general parameters ----------------------------------------------------------

set.seed(123)

# simulation time steps
sim_days <- 100 

#################################################
# life history parameters
#################################################

# Number of patches
patches = 5     

# Initial number of individuals per patch
initial_pop = 1000                  # carrying capacity
n_per_patch <- create_n_per_patch(patches, initial_pop)   # Initial number of individuals per patch
sex_alleles <- c("X", "Y") 
prob_wildtype2 = 0.05 
prob_wildtype1 = 1 - prob_wildtype2 
stages <- c("egg", "larva", "pupa", "adult")
dd_effect = 0.00001
max_survival = 0.85


# the adult male population size at which the daily probability of mating is 
# 0.5 (North and Godfray, Malar J (2018) 17:140)
beta <- 200      

# Probability that a female find a blood meal
bloodmeal_prob <- 0.40                   

# low degree-day threshold
ldt <- c(egg = 8.19, 
         larva = 9.33, 
         pupa = 10.3)

# dispersal parameters

# lambda controls the rates at which probability between patches decreases with distances 
lambda <- 0.1

dispersal_prop <- 0.0005


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
# Gene drive parameters
#################################################


n_loci <- 10   # number of loci
release_freq = 0.1     # drive release frequency: percentage of adults
cut_rate = 0.98         # gene cut_rate probability 
homing_rate =  0.98
decay <- 0.5            #controls the rate at which the covariance between two loci decreases with distance


################################################
# environmental parameters
################################################

# create coordinates and dipersal matrix for the patches/locations 

# random location coordinates
set.seed(2025)

# coords <- as.data.frame(100 * matrix(runif(patches * 2), ncol = 2))
# colnames(coords) <- c("x","y")

coords <- data.frame(
  x = seq(0, 100, length.out = patches),
  y = rep(50, patches)
)


plot(coords, cex = 10)
text(coords, labels = 1:patches)


# dispersal matrix
neg_exponet_model <- metapop (coords = coords, 
                                 lambda = lambda, 
                                 disp_prob = dispersal_prop)

adjacency_matrix <- step_stone(n_patches = patches, 
                               disp_prob = dispersal_prop)




# This bit of code generates random daily temperature and humidity to estimates  
# growth degree-day required for each stage during transition and survival

# temp_max  <-  matrix(runif(patches * sim_days, 
#                         min = 30, max = 38), 
#                   nrow = sim_days, ncol = patches) 
# 
# temp_min <-  matrix(runif(patches * sim_days, 
#                        min = 16, max = 22), 
#                  nrow = sim_days, ncol = patches)
# 
# humidity <- matrix(rtruncnorm(patches * sim_days, a = 0, b = 100,
#                         mean = 80, sd = 10), 
#                   nrow = sim_days, ncol = patches)

# Or as suggested to prevent the significant impact on the population dynamics,
# use a fixed value instead 

temp_max <- 30
temp_min <- 25
humidity <- 83


# size of the habitat in cm^2 to estimate the size of the habitat for aquatic 
# stages (check with Nick)
s_area <- 3333.2                   





