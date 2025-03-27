library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse)

##TO DO!
# introduce a functional exponential dispersal model. DONE!
# implement density-dependence survival at the larval stage. DONE!
# introduce effect of temperature and humidity on transition/growth and survival
# introduce gene drive model 

###########################################
#               PARAMETERS                #
###########################################

set.seed(21032025)

patches <- 5                        # Number of patches
n_per_patch <- c(500, 
                 150, 
                 10, 
                 5, 
                 5)                 # Initial number of individuals per patch

daily_survival <- c(egg = 0.8,      # daily survival prob
                    larva = 0.8, 
                    pupa = 0.8, 
                    adult = 0.8)    
daily_transition <- c(egg = 0.55,
                      larva = 0.6,
                      pupa = 0.7)   # daily transition prob
growth_rate <- 5                    # Number of offspring per day per female mosquitoe 
carry_k <- 2000                   # Carrying capacity
mate_prob <- 0.65                   # Probability of mating
sim_days <- 20                      # Number of simulation in days
bloodmeal_prob <- 0.65              # Probability that a female find a blood meal
dispersal_frac <- 0.2



###########################################
#  INITIALISE POPULATION WITH ATTRIBUTES  #
###########################################


ini_pop <- function(patches, n_per_patch) {
  patches_pop <- list()
  
  for (i in 1:patches) {
      x_value <- runif(1, 0, 100)
      y_value <- runif(1, 0, 100)
      patches_pop[[i]] <- tibble(
      patch = i,
      sex = rbinom(n_per_patch[i], 1, 0.5), # Female == 1, random sex
      stage = sample(c("egg", "larva", "pupa", "adult"), n_per_patch[i], replace = TRUE),
      alive = TRUE,
      x = rep(x_value, n_per_patch[i]),
      y = rep(y_value, n_per_patch[i])
     )
  }
  
  return(patches_pop)
}


#Check
pop <- ini_pop(patches, n_per_patch)

sum((sapply(pop, nrow)))


############################################
#      GROWTH (including reproduction)     #
############################################


growth <- function(pop_patches, 
                   mate_prob, 
                   bloodmeal_prob, 
                   growth_rate,
                   daily_survival,
                   daily_transition,
                   carry_k) { 
  
  updated_pop_patches <- list()
  
  for (i in seq_along(pop_patches)) {
    pop <- pop_patches[[i]]  
    
    n_fem <- pop |> filter(sex == 1, stage == "adult")       # All females
    mated_fem <- rbinom(1, nrow(n_fem), mate_prob)           # All mated females
    bloodfed_fem <- rbinom(1, mated_fem, bloodmeal_prob)     # Probability that a female finds a blood meal

    exp_offspring <- bloodfed_fem * growth_rate
    n_offspring <- rpois(1, exp_offspring)  
    
    offspring <- tibble()
    
    if (n_offspring > 0) {
      
        offspring <- tibble(
        patch = pop$patch[1],
        sex = rbinom(1, n_offspring, 0.5),  
        stage = rep("egg", n_offspring),
        alive = TRUE,
        x = pop$x[1],
        y = pop$y[1]
     )
    
    
    pop <- bind_rows(pop, offspring)
      
    }
   
    # Density-dependent survival for larval stage
    larva_count <- sum(pop$stage == "larva")
    density_dependent_survival <- daily_survival["larva"] * (1 - (larva_count / carry_k))
    
    # Probability of survival in each time step  
    pop <- pop |>
      mutate(
        alive = case_when(
          stage == "egg" ~ rbinom(n(), 1, daily_survival["egg"]),
          stage == "larva" ~ rbinom(n(), 1, density_dependent_survival),
          stage == "pupa" ~ rbinom(n(), 1, daily_survival["pupa"]),
          stage == "adult" ~ rbinom(n(), 1, daily_survival["adult"]),
          TRUE ~ NA_integer_
        ),
        alive = alive == 1
      )   
   
     # stage transition probability per time step 
      pop <- pop |>
      mutate(
        stage = case_when(
        stage == "egg" & rbinom(n(), 1, daily_transition["egg"]) == 1 ~ "larva",
        stage == "larva" & rbinom(n(), 1, daily_transition["larva"]) == 1 ~ "pupa",
        stage == "pupa" & rbinom(n(), 1, daily_transition["pupa"]) == 1 ~ "adult",
        .default = stage)
      )
    
    pop <- filter(pop, alive)
      
    updated_pop_patches[[i]] <- pop
  }
  return(updated_pop_patches)
}
  

# # check
grown_pop <- growth(pop_patches = pop,
                    mate_prob = 0.6,
                   bloodmeal_prob = 0.7,
                   growth_rate = 5,
                   daily_survival = daily_survival,
                   daily_transition = daily_transition,
                   carry_k = 2000)

# sum((sapply(grown_pop, nrow)))



###########################################
#                DISPERSAL                #
###########################################


dispersal <- function(pop, dispersal_frac) {
  
  # build a function for distance matrix based of on the population x & y coords.
  
  create_dist_matrix <- function(pop) {
    # This ensures we get the coordinates of the patches
    coords <- pop |>
      bind_rows() |>
      select(patch, x, y) |>
      distinct()            
    # Calculate pairwise Euclidean distances between patches
    dist_matrix <- as.matrix(dist(coords[, c("x", "y")], method = "euclidean"))
    # Return the distance matrix
    return(dist_matrix)
  }
  
  
  #create distance matrix using function
  dist_matrix <- create_dist_matrix(pop)
  
  # create exponential dispersal kernel using the distance matrix. 
  lambda <- 1/mean(dist_matrix) # range is = 1/average distance in the matrix
  
  dispersal_kernel <- exp(-lambda * dist_matrix)
  
  # set the diagonal elements to 0 to prevent self-dispersal
  
  diag(dispersal_kernel) <- 0
  
  
  # make these columns sum to 1 to get probability of moving to other patch
  # *if* they left. This dispersal matrix gives the probability of the vector
  # vector moving between patches
  rel_dispersal_matrix <- sweep(dispersal_kernel, 1,
                                rowSums(dispersal_kernel), FUN = "/")
  
  # sum(rel_dispersal_matrix[2,])
  
  
  # normalise these to have the overall probability of dispersing to that patch,
  # and add back the probability of remaining
  dispersal_matrix <- dispersal_frac * rel_dispersal_matrix +
    (1 - dispersal_frac) * diag(nrow(dispersal_kernel))  
  
  dispersed_pop <- pop

  for (i in 1:length(pop)) {
    patch <- pop[[i]]

    # adults in the patch capable of dispersing
    adults <- patch[patch$stage == "adult" & patch$alive == TRUE, ]

    # If there are no adults, skip dispersal for this patch
    if (nrow(adults) == 0) next
    
    # Get the dispersal probabilities for this individual according to the dispersal matrix
    dispersal_probs <- dispersal_matrix[i,]

    for (j in 1:nrow(adults)) {

      individual <- adults[j, ]
      
      # Use multinomial distribution to sample a new patch
      new_patch_index <- which(rmultinom(1, 1, dispersal_probs) == 1)
      
      # Add the individual to the selected new patch
      dispersed_pop[[new_patch_index]] <- rbind(dispersed_pop[[new_patch_index]], individual)
      
      # Remove the individual from the original patch
      dispersed_pop[[i]] <- dispersed_pop[[i]][!rownames(dispersed_pop[[i]]) %in% rownames(adults)[j], ]
    }
  }
  return(dispersed_pop)
}



# # check
# disp_pop <- dispersal(grown_pop, dispersal_frac)
# sum((sapply(disp_pop, nrow)))


###########################################
#               SIMULATION                #
###########################################


simulation <- function(patches, 
                       n_per_patch, 
                       mate_prob, 
                       bloodmeal_prob, 
                       growth_rate, 
                       daily_survival, 
                       daily_transition,
                       sim_days,
                       carry_k,
                       dispersal_frac) {
  
  pop <- ini_pop(patches, n_per_patch)
  
  patch_sizes <- list()
  stage_distributions <- list()

  for (day in 1:sim_days) {
    cat("Day", day, "Completed\n")
    
    
    # Growth with reproduction
    pop <- growth(pop,
                  mate_prob,
                  bloodmeal_prob,
                  growth_rate,
                  daily_survival,
                  daily_transition,
                  carry_k)

    # Dispersal
    pop <- dispersal(pop, dispersal_frac)
    
    # Track daily population sizes per patch
    patch_sizes[[day]] <- sapply(pop, nrow)  

    # Track stage distribution per patch
    stage_distributions[[day]] <- lapply(pop, function(patch_pop) {
      table(patch_pop$stage)  
    })
  }
  
  # Return the collected data
  list(
    patch_sizes = patch_sizes,
    stage_distributions = stage_distributions,
    final_pop = pop  
  )
}


###########################################
#            RUN   SIMULATION             #
###########################################

sim <- simulation(patches = 5, 
                  n_per_patch = n_per_patch, 
                  mate_prob = 0.6, 
                  bloodmeal_prob = bloodmeal_prob, 
                  growth_rate = growth_rate, 
                  daily_survival = daily_survival, 
                  daily_transition = daily_transition,
                  carry_k = 2000,
                  sim_days = 20,
                  dispersal_frac= 0.2)



