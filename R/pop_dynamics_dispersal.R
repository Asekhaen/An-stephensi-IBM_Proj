#########
#simple population dynamics and dispersal model ####


library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr)

# Assign (random0 coordinates to the locations 
n_locs = 15
coords <- 100 * matrix(runif(n_locs * 2), ncol = 2)
k_capacity <- 10000
growth_rate <- 4


#create population matrix and initialise vector population 

initilise_pop <- function(n_locs, 
                          coords,
                          k_capacity){
 
  # invasive status of each location: 1 = native, 0 = invasive
  status <- sample(c(1,0), n_locs, replace = TRUE)   
  
  # vector population: populate each location based on their invasion status by  
  # multiplying each location "status" with the location carrying capacity which 
  # is currently fixed.
  vector_pop <- status * k_capacity

  # Bind parameters together
  pop <- cbind(coords, status, vector_pop)
  
  # Assign column names
  colnames(pop) <- c("x", "y", "status", "vector_pop")
  
  return(pop)
}



pop <-initilise_pop(n_locs = n_locs, coords = coords, k_capacity = k_capacity)




#### Beverton_Holt population growth with density dependence #####

dd_growth <- function(pop, growth_rate, k_capacity) {
  for (i in 1:nrow(pop)) {
    # Only update population if the population is greater than 0
    if (pop[i, "vector_pop"] > 0) {
      # Beverton-Holt growth model
      pop[i, "vector_pop"] <- round(pop[i, "vector_pop"] * growth_rate / 
                                      (1 + (pop[i, "vector_pop"] * 
                                              (growth_rate - 1)) / k_capacity))
    }
  }
  return(pop)
}


grown_state <- dd_growth(pop = pop, growth_rate = growth_rate, k_capacity)

##### dispersal function ####
# the dispersal function in the  DSDM model consist of three dispersal approahces
# to capture all possible dispersal routes including long and short distance
# dispersal and those resulting from human movement etc. The current model
# implememted exponential dispersal only.

# build a euclidean matrix (distance matrix)) based on the location coordinates
distance_matrix <- as.matrix(dist(coords))


# model exponential dispersal kernel using the distance matrix. This gives the
# probability of a vector dispersal between a pair of location
lambda <- 1/mean(distance_matrix)

exponential_dispersal_kernel <- exp(-lambda * distance_matrix)

# set the diagonal elements to 0
diagonal <- diag(nrow = n_locs, ncol = n_locs)
exponential_dispersal_kernel <- exponential_dispersal_kernel * (1 - diagonal)

# make these columns sum to 1 to get probability of moving to each other patch
# *if* they left. This dispersal matrix gives the probability of the vector
# vector moving between patches
rel_dispersal_matrix <- sweep(exponential_dispersal_kernel, 1,
                              rowSums(exponential_dispersal_kernel), FUN = "/")

#sum(rel_dispersal_matrix[1,])

# probability of dispersion. This gives the probability of the fraction of the
#vectors leaving a patch/location. (fixed for now, but can be modify later based
#on knowledge of the dispersal probability and vector species)
disp_fraction <- 0.02

# normalise these to have the overall probability of dispersing to that patch,
# and add back the probability of remaining
dispersal_matrix <- disp_fraction * rel_dispersal_matrix +
  (1 - disp_fraction) * diagonal



# Dispersal function, updates population based on dispersal probabilities (Stochastic1)

dispersal <- function(pop, dispersal_matrix) {
  # Loop over each location
  for (i in 1:nrow(pop)) {
    # For each location, calculate the number of individuals dispersing to other loc.
    for (j in 1:nrow(pop)) {
      if (i != j) {  # Only disperse to other patches, not itself
        # Calculate the number of individuals to disperse from patch i to patch j
        dispersing_vector <- rbinom(1, pop[i, "vector_pop"], dispersal_matrix[i, j])
        
        # Update populations: subtract from source patch, add to target patch
        pop[i, "vector_pop"] <- pop[i, "vector_pop"] - dispersing_vector
        pop[j, "vector_pop"] <- pop[j, "vector_pop"] + dispersing_vector
      }
    }
  }
  return(pop)
}

disperesed_pop <- dispersal(pop = grown_state, dispersal_matrix = dispersal_matrix)

# Dispersal function, updates population based on dispersal probabilities (Deterministic)

# dispersal <- function(pop, dispersal_matrix) {
#   # Loop over each patch (location)
#   for (i in 1:nrow(pop)) {
#     # For each patch, calculate the number of individuals dispersing to other patches
#     for (j in 1:nrow(pop)) {
#       if (i != j) {  # Only disperse to other patches, not itself
#         # Calculate the number of individuals to disperse from patch i to patch j
#         dispersing_vector <- round(pop[i, "vector_pop"] * dispersal_matrix[i, j])
# 
#         # Update populations: subtract from source patch, add to target patch
#         pop[i, "vector_pop"] <- pop[i, "vector_pop"] - dispersing_vector
#         pop[j, "vector_pop"] <- pop[j, "vector_pop"] + dispersing_vector
#       }
#     }
#   }
#   return(pop)
# }





# Simulation Function - Corrected to store population at each location and time step
sim <- function(n_locs, coords, k_capacity, growth_rate, time_step, dispersal_matrix) {
  
  # Store population at each location for each time step
  pop <- initilise_pop(n_locs, coords, k_capacity)
  
  # Matrix to store population at each time step
  pop_per_loc <- matrix(NA, nrow = time_step, ncol = n_locs)
  
  # Loop over time steps
  for (t in 1:time_step) {
    cat("Time", t, "\n")
    
    # Apply Beverton-Holt growth model
    pop <- dd_growth(pop, growth_rate, k_capacity)
    
    # Apply dispersal
    pop <- dispersal(pop, dispersal_matrix)
    
    # Store population at each location at the current time step
    pop_per_loc[t, ] <- pop[, "vector_pop"]
  }
  
  return(pop_per_loc)
}





pop_sim <- sim(n_locs = n_locs, 
               coords = coords,
               k_capacity = 5000, 
               growth_rate = growth_rate,
               time_step = 20,
               dispersal_matrix = dispersal_matrix)
                   



# Convert the returned population matrix into a dataframe for plotting

pop_df <- as.data.frame(pop_sim)
pop_df$time_step <- 1:nrow(pop_df) # Add the time steps as a column

# Reshape the data to long format for ggplot2
pop_long <- pop_df |>
  pivot_longer(cols = -time_step, names_to = "location", values_to = "population_size")

# Plot the population size over time for each location
ggplot(pop_long, aes(x = time_step, y = population_size, color = location)) +
  geom_line() +
  labs(title = "Population dynamics",
       x = "Time Step",
       y = "Population Size",
       color = "patch") +
  theme_minimal() +
  theme(legend.position = "right")






