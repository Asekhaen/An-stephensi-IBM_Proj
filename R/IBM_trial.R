library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)



###########################################
#               PARAMETERS                #
###########################################

patches <- 3                        # Number of patches
num_per_patch <- c(100,10,10)       # Initial number of individuals per patch
num_loci <- 4                       # Number of loci
dispersal_rate <- 0.01              # Base dispersal rate (adjust this if necessary)
daily_survival <- c(egg = 0.8,      # daily survival prob
                    larva = 0.8, 
                    pupa = 0.8, 
                    adult = 0.8)    
daily_transition <- c(egg = 0.55,
                      larva = 0.6,
                      pupa = 0.7)   # daily transition prob
offspring_day <- 5                  # Number of offspring per day per female mosquitoes 
carry_k <- 10000                    # Carrying capacity
mate_prob <- 0.65                   # Probability of mating
sim_days <-150                      # Number of simulation in days
bloodmeal_prob <- 0.65              # Probability that a female find a blood meal



###########################################
#  INITIALISE POPULATION WITH ATTRIBUTES  #
###########################################


# Function to create a chromosome
create_chromosome <- function(num_loci) {
  rbinom(num_loci, 1, 0.9) # For example; Drive = 1, wild = 0
}


ini_pop <- function(patches, num_per_patch, num_loci) {
  patches_pop <- list()
  
  for (i in 1:patches) {
    patches_pop[[i]] <- tibble(
      age = 0,                           # Age of each individual in the patch
      sex = rbinom(num_per_patch[i], 1, 0.5), # Female == 1, random sex
      stage = "egg",
      alive = TRUE,
      chromosome1 = lapply(1:num_per_patch[i], function(x) create_chromosome(num_loci)), # Chromosome 1
      chromosome2 = lapply(1:num_per_patch[i], function(x) create_chromosome(num_loci))  # Chromosome 2
    )
  }
  
  return(patches_pop)
}





###########################################
#               REPRODUCTION              #
###########################################


#Beverton-Holt density-dependent 

bev_holt <- function(N, offspring_day, carry_k) {
  return((offspring_day * N) / (1 + (offspring_day - 1) * N / carry_k))
}



reprod <- function(pop_patches, mate_prob, bloodmeal_prob, offspring_day, carry_k) { 
  
  updated_pop_patches <- list()
  
  for (i in seq_along(pop_patches)) {
    pop <- pop_patches[[i]]  
    
    num_fem <- sum(pop$sex == 1 & pop$stage == "adult")  # All females
    num_mated_fem <- sum(runif(num_fem) < mate_prob)  # All mated females
    bloodfed_fem <- sum(runif(num_mated_fem) < bloodmeal_prob) #probability that a female finds a blood meal
    num_males <- sum(pop$sex == 0 & pop$stage == "adult")  # All males
    
    n_pairs <- min(bloodfed_fem, num_males)
    
    expected_off <- bev_holt(n_pairs, offspring_day, carry_k)  
    actual_off <- rpois(1, expected_off)  
    
    offspring <- tibble()
    
    if (actual_off > 0 && n_pairs > 0) {
      
      f_parent <- pop |>
        filter(sex == 1 & stage == "adult") |>
        sample_n(n_pairs, replace = FALSE)
      
      m_parent <- pop |>
        filter(sex == 0 & stage == "adult") |>
        sample_n(n_pairs, replace = FALSE)
      
      pairs <- rep(1:n_pairs, length.out = actual_off)
      
      offspring <- tibble(
        age = rep(0, actual_off),
        sex = rbinom(actual_off, 1, 0.5),  
        stage = rep("egg", actual_off),
        alive = TRUE,
        
        chromosome1 = ifelse(runif(actual_off) < 0.5, f_parent$chromosome1[pairs], 
                             m_parent$chromosome1[pairs]),
        chromosome2 = ifelse(runif(actual_off) < 0.5, f_parent$chromosome2[pairs], 
                             m_parent$chromosome2[pairs])
     )
    }
    
    pop <- bind_rows(pop, offspring)
    
    updated_pop_patches[[i]] <- pop
  }
  
  return(updated_pop_patches)
}





###########################################
#                 GROWTH                  #
###########################################



growth <- function(updated_pop, daily_survival, daily_transition){
  aged_pop <- list()
  
  for (i in seq_along(updated_pop)) {
    pop <- updated_pop[[i]] 
    
#daily survival of vector stages 
    pop <- pop |>
      mutate(
        alive = case_when(
          stage == "egg" ~ rbinom(n(), 1, daily_survival["egg"]),
          stage == "larva" ~ rbinom(n(), 1, daily_survival["larva"]),
          stage == "pupa" ~ rbinom(n(), 1, daily_survival["pupa"]),
          stage == "adult" ~ rbinom(n(), 1, daily_survival["adult"]),
        ),
        alive = alive == 1
      ) |>
      filter(alive)
    
    
#Density-dependence growth for larval stage
    
    
  
#daily transition probability between stages
  pop <- pop |>
      mutate(stage = case_when(
        stage == "egg" & rbinom(n(), 1, daily_transition["egg"]) == 1 ~ "larva",
        stage == "larva" & rbinom(n(), 1, daily_transition["larva"]) == 1 ~ "pupa",
        stage == "pupa" & rbinom(n(), 1, daily_transition["pupa"]) == 1 ~ "adult",
        .default = stage
      ),
      age = age +1)

      pop <- pop |> filter(alive)
   
       aged_pop[[i]] <- pop
  }
  
  return(aged_pop)
}





###########################################
#                DISPERSAL                #
###########################################


# Generate dispersal matrix based on the distance matrix
dispersal_matrix <- dispersal_kernel(distance_matrix)

# Ensure no individual disperses back to the same patch
diag(dispersal_matrix) <- 0  # Make sure diagonal elements are 0 to prevent self-dispersal

# Function for dispersal of individuals
Dispersal <- function(aged_pop, dispersal_matrix, dispersal_rate) {
  dispersed_pop <- list()
  # Loop over each patch
  for (i in 1:length(aged_pop)) {
    patch <- aged_pop[[i]]
    
    # Identify adults in the patch
    adults <- patch |> filter(stage == "adult" & alive == TRUE & runif(n()) < dispersal_rate)
    
    # If there are no adults, skip dispersal for this patch
    if (nrow(adults) == 0) next
    
    # Loop over each adult individual in the patch
    for (j in 1:nrow(adults)) {
      # Get the dispersal probabilities for this individual (based on the distance to other patches)
      dispersal_probs <- dispersal_matrix[i,]
      
      # Normalize the dispersal probabilities to make sure they sum to 1
      dispersal_probs <- dispersal_probs / sum(dispersal_probs)
      
      # Simulate dispersal (choose a patch based on the dispersal probabilities)
      new_patch <- sample(1:length(aged_pop), size = 1, prob = dispersal_probs)
      
      # Move the individual to the new patch (we're keeping track of the movement, not modifying the actual data here)
      aged_pop[[new_patch]] <- bind_rows(aged_pop[[new_patch]], adults[j,])
      
      # Remove the individual from the current patch (if they moved)
      aged_pop[[i]] <- aged_pop[[i]] |>
        filter(row_number() != j)
    }
    dispersed_pop[[i]] <- aged_pop  
  }
  return(aged_pop)
}




###########################################
#              SIMULATION                 #
###########################################


simulation <- function(sim_days, patches, num_per_patch, max_age, 
                       num_loci, mate_prob, bloodmeal_prob, offspring_day, carry_k, 
                       daily_survival, daily_transition, dispersal_rate) {
  
  pop <- ini_pop(patches, num_per_patch, num_loci)
  
  patch_sizes <- list()
  age_distributions <- list()
  stage_distributions <- list()
  sex_distributions <- list()
  
  for (day in 1:sim_days) {
    cat("Day", day, "\n")
    
    # Reproduction
    updated_pop <- reprod(pop, mate_prob, bloodmeal_prob, offspring_day, carry_k)
    
    # Growth
    aged_pop <- growth(updated_pop, daily_survival, daily_transition)
    
    # Dispersal
    dispersed_pop <- Dispersal(aged_pop, dispersal_matrix, dispersal_rate)
    
    
    pop <- dispersed_pop   # Update population for the next day
    
    # Track daily population sizes per patch
    daily_sizes <- sapply(pop, nrow)  
    patch_sizes[[day]] <- daily_sizes
    
    # Track age distribution per patch
    daily_age_dist <- lapply(pop, function(patch_pop) {
      table(patch_pop$age)  
    })
    age_distributions[[day]] <- daily_age_dist
    
    #Track sex distribution 
    
    daily_sex_dist <- lapply(pop, function(patch_pop) {
      table(patch_pop$sex)
    })
    sex_distributions[[day]] <- daily_sex_dist
    
    
    # Track stage distribution per patch
    daily_stage_dist <- lapply(pop, function(patch_pop) {
      table(patch_pop$stage)  
    })
    stage_distributions[[day]] <- daily_stage_dist
  }
  
  # Return the collected data
  list(
    patch_sizes = patch_sizes,
    age_distributions = age_distributions,
    stage_distributions = stage_distributions,
    sex_distributions = sex_distributions,
    final_pop = pop  
  )
}




###########################################
#            RUN   SIMULATION             #
###########################################

sim <- simulation(sim_days, 
                  patches, 
                  num_per_patch, 
                  max_age, 
                  num_loci, 
                  mate_prob, 
                  bloodmeal_prob, 
                  offspring_day, 
                  carry_k, 
                  daily_survival, 
                  daily_transition, 
                  dispersal_rate)




