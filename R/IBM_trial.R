library(tibble)
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

set.seed(04042025)


source("R/parameters.R")

###########################################
#  INITIALISE POPULATION WITH ATTRIBUTES  #
###########################################

# create corrdinates for the patches/locations 
coords <- as.data.frame(100 * matrix(runif(patches * 2), ncol = 2))
colnames(coords) <- c("x","y")


# initialise population
ini_pop <- function(patches, n_per_patch, coords) {
  patches_pop <- list()
  
  for (i in 1:patches) {
      patches_pop[[i]] <- tibble(
      patch = i,
      sex = rbinom(n_per_patch[i], 1, 0.5), # Female == 1, random sex
      stage = sample(c("egg", "larva", "pupa", "adult"), n_per_patch[i], replace = TRUE),
      alive = TRUE,
      x = rep(coords$x[i], n_per_patch[i]),
      y = rep(coords$y[i], n_per_patch[i])
     )
  }
  
  return(patches_pop)
}


#Check
pop <- ini_pop(patches, n_per_patch, coords)

sum((sapply(pop, nrow)))


############################################
#      GROWTH (including reproduction)     #
############################################


growth <- function(pop_patches, 
                   mate_prob, 
                   bloodmeal_prob, 
                   fecundity,
                   daily_survival,
                   daily_transition,
                   carry_k) { 
  
  updated_pop_patches <- list()
  
  for (i in seq_along(pop_patches)) {
    pop <- pop_patches[[i]]  
    
    n_fem <- pop |> filter(sex == 1, stage == "adult")       # All females
    mated_fem <- round(rbinom(1, nrow(n_fem), mate_prob))           # All mated females
    bloodfed_fem <- round(rbinom(1, mated_fem, bloodmeal_prob))     # Probability that a female finds a blood meal

    exp_offspring <- bloodfed_fem * fecundity
    n_offspring <- rpois(1, exp_offspring)  
    
    offspring <- tibble()
    
    if (n_offspring > 0) {
      
        offspring <- tibble(
        patch = pop$patch[1],
        sex = rbinom(n_offspring,1, 0.5),  
        stage = rep("egg", n_offspring),
        alive = TRUE,
        x = pop$x[1],
        y = pop$y[1]
     )
    
    
    pop <- bind_rows(pop, offspring)
      
    }
   
  # Density-dependent survival for larval stage
    larva_count <- sum(pop$stage == "larva")
    # density_dependent_survival <- daily_survival["larva"] * (1 - (larva_count / carry_k))
    
    density_dependence <- 1 - larva_count / carry_k
    density_dependence <- max(density_dependence, 0)
    density_dependent_survival <- daily_survival["larva"] * density_dependence
    
    # Probability of survival in each time step (this approach produced NAs 
    # because of patches with empty stages) 
    # pop <- pop |> mutate(
    #     alive = case_when(
    #       stage == "egg" ~ rbinom(n(), 1, daily_survival["egg"]),
    #       stage == "larva" ~ rbinom(n(), 1, density_dependent_survival),
    #       stage == "pupa" ~ rbinom(n(), 1, daily_survival["pupa"]),
    #       stage == "adult" ~ rbinom(n(), 1, daily_survival["adult"]),
    #       TRUE ~ NA_integer_
    #     ),
    #     alive = alive == 1
    #   )
   
    pop <- pop |>
      mutate(alive = case_when(
        stage == "egg" ~ runif(n()) < daily_survival["egg"],
        stage == "larva" ~ runif(n()) < density_dependent_survival,
        stage == "pupa" ~ runif(n()) < daily_survival["pupa"],
        stage == "adult" ~ runif(n()) < daily_survival["adult"]
      ) & alive)
    
    
     # stage transition probability per time step 
      pop <- pop |>
      mutate(
        stage = case_when(
        stage == "egg" & rbinom(n(), 1, daily_transition["egg"]) == 1 ~ "larva",
        stage == "larva" & rbinom(n(), 1, daily_transition["larva"]) == 1 ~ "pupa",
        stage == "pupa" & rbinom(n(), 1, daily_transition["pupa"]) == 1 ~ "adult",
        TRUE ~ stage
        )
      )
    
    pop <- filter(pop, alive)
      
    updated_pop_patches[[i]] <- pop
  }
  return(updated_pop_patches)
}
  

# # check
grown_pop <- growth(pop_patches = pop,
                    mate_prob = 0.81,
                   bloodmeal_prob = 0.7,
                   fecundity = 9,
                   daily_survival = daily_survival,
                   daily_transition = daily_transition,
                   carry_k = 1000)

# sum((sapply(grown_pop, nrow)))



###########################################
#                DISPERSAL                #
###########################################

# build a function for distance matrix based of on the population x & y coords.

source("R/dispersal_matrix.R") 

dispersal <- function(pop, dispersal_matrix) {

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

      # Update the coordinates of the individual to match the new patch's coordinates
      dispersed_pop[[new_patch_index]]$x[nrow(dispersed_pop[[new_patch_index]])] <- coords$x[new_patch_index]
      dispersed_pop[[new_patch_index]]$y[nrow(dispersed_pop[[new_patch_index]])] <- coords$y[new_patch_index]
      
      # Update the patch number for the dispersed individual
      dispersed_pop[[new_patch_index]]$patch[nrow(dispersed_pop[[new_patch_index]])] <- new_patch_index

      # Remove the individual from the original patch
      dispersed_pop[[i]] <- dispersed_pop[[i]][!rownames(dispersed_pop[[i]]) %in% rownames(adults)[j], ]
    }
  }
  return(dispersed_pop)
}



# # check
disp_pop <- dispersal(pop = grown_pop, dispersal_matrix)
sum((sapply(disp_pop, nrow)))



###########################################
#               SIMULATION                #
###########################################


simulation <- function(patches,
                       coords = coords,
                       n_per_patch, 
                       mate_prob, 
                       bloodmeal_prob, 
                       fecundity, 
                       daily_survival, 
                       daily_transition,
                       sim_days,
                       carry_k,
                       dispersal_matrix) {
  
  pop <- ini_pop(patches, n_per_patch, coords)
  
  patch_sizes <- list()
  stage_distributions <- list()

  for (day in 1:sim_days) {
    cat("Day", day, "Completed\n")
    
    
    # Growth with reproduction
    pop <- growth(pop,
                  mate_prob,
                  bloodmeal_prob,
                  fecundity,
                  daily_survival,
                  daily_transition,
                  carry_k)

    # Dispersal
    pop <- dispersal(pop, dispersal_matrix)
    
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

sim <- simulation(patches = 6, 
                  coords = coords,
                  n_per_patch = n_per_patch, 
                  mate_prob = 0.81, 
                  bloodmeal_prob = 0.7, 
                  fecundity = 9, 
                  daily_survival = daily_survival, 
                  daily_transition = daily_transition,
                  carry_k = 1000,
                  sim_days = 20,
                  dispersal_matrix = dispersal_matrix)
