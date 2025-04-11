# Functions for running the stephensii IBM

# initialise population
ini_pop <- function(patches, n_per_patch, coords, loci) {
  patches_pop <- list()
  
  for (i in 1:patches) {
    patches_pop[[i]] <- tibble(
      patch = i,
      sex = rbinom(n_per_patch[i], 1, 0.5), # Female == 1, random sex
      stage = sample(c("egg", "larva", "pupa", "adult"), n_per_patch[i], replace = TRUE),
      alive = TRUE,
      x = rep(coords$x[i], n_per_patch[i]),
      y = rep(coords$y[i], n_per_patch[i]),
      allele1 = matrix(sample(c(0, 1), n_per_patch[i] * n_loci, replace = TRUE, prob = c(0.8,0.2)), ncol = n_loci,), # 0 = wild-type, 1 = drive allele
      allele2 = matrix(sample(c(0, 1), n_per_patch[i] * n_loci, replace = TRUE, prob = c(1,0)), ncol = n_loci)
    )
  }
  
  return(patches_pop)
}

growth <- function(pop_patches, 
                   mate_prob, 
                   bloodmeal_prob, 
                   fecundity,
                   drive_conversion_prob,
                   n_loci,
                   daily_survival,
                   daily_transition,
                   carry_k) { 
  updated_pop_patches <- list()
  
  for (i in seq_along(pop_patches)) {
    pop <- pop_patches[[i]]  
    
    n_male <- pop |> filter(sex == 0, stage == "adult")       # All females
    n_fem <- pop |> filter(sex == 1, stage == "adult")            # All females
    mated_fem <- round(rbinom(1, nrow(n_fem), mate_prob))         # All mated females
    bloodfed_fem <- round(rbinom(1, mated_fem, bloodmeal_prob))   # Probability that a female finds a blood meal
    
    exp_offspring <- bloodfed_fem * fecundity
    n_offspring <- rpois(1, exp_offspring)  
    
    offspring <- tibble()
    
    
    # # Genetic inheritance and drive conversion with no effect
    # if (n_offspring > 0) {
    # 
    # offspring <- tibble(
    #   patch = pop$patch[1],
    #   sex = rbinom(n_offspring,1, 0.5),
    #   stage = rep("egg", n_offspring),
    #   alive = TRUE,
    #   x = pop$x[1],
    #   y = pop$y[1],
    #   allele1 = matrix(sample(c(n_fem$allele1, n_fem$allele2), n_offspring * n_loci, replace = TRUE), ncol = n_loci),
    #   allele2 = matrix(sample(c(n_male$allele1, n_male$allele2), n_offspring * n_loci, replace = TRUE), ncol = n_loci)
    # )
    
    
    # Genetic inheritance and drive conversion
    # if (n_offspring > 0) { 
    if (n_offspring > 0 && nrow(n_male$allele1) > 0 && nrow(n_male$allele2) > 0) {
      allele1_offspring <- matrix(sample(c(n_fem$allele1, n_fem$allele2), n_offspring * n_loci, replace = TRUE), ncol = n_loci)
      allele2_offspring <- matrix(sample(c(n_male$allele1, n_male$allele2), n_offspring * n_loci, replace = TRUE), ncol = n_loci)
      
      for (j in 1:n_offspring) {
        for (k in 1:n_loci) {
          if (allele1_offspring[j, k] == 1) {
            allele2_offspring[j, k] <- ifelse(runif(1) < drive_conversion_prob, 1, allele2_offspring[j, k])
          }
          if (allele2_offspring[j, k] == 1) {
            allele1_offspring[j, k] <- ifelse(runif(1) < drive_conversion_prob, 1, allele1_offspring[j, k])
          }
        }
      }
      
      offspring <- tibble(
        patch = pop$patch[1],
        sex = rbinom(n_offspring, 1, 0.5),
        stage = rep("egg", n_offspring),
        alive = TRUE,
        x = pop$x[1],
        y = pop$y[1],
        allele1 = allele1_offspring,
        allele2 = allele2_offspring
      )
      
      # Gene drive lethal effect: 100% mortality of individual with homozygous loci
      
      for (j in 1:nrow(offspring)) {
        for (k in 1:n_loci) {
          if (offspring$allele1[j, k] == 1 && offspring$allele2[j, k] == 1) {
            offspring$alive[j] <- FALSE
          }
        }
      }
      
      pop <- bind_rows(pop, offspring)
      
    }
    
    
    # Density-dependent survival for larval stage
    larva_count <- sum(pop$stage == "larva")
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

#  Dispersal function 
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
      dispersed_pop[[new_patch_index]] <- bind_rows(dispersed_pop[[new_patch_index]], individual) 
      
      # Update the coordinates of the individual to match the new patch's coordinates
      dispersed_pop[[new_patch_index]]$x[nrow(dispersed_pop[[new_patch_index]])] <- coords$x[new_patch_index]
      dispersed_pop[[new_patch_index]]$y[nrow(dispersed_pop[[new_patch_index]])] <- coords$y[new_patch_index]
      
      # Update the patch number for the dispersed individual
      dispersed_pop[[new_patch_index]]$patch[nrow(dispersed_pop[[new_patch_index]])] <- new_patch_index
      
      # Remove the individual from the original patch
      dispersed_pop[[i]] <- dispersed_pop[[i]][-which(rownames(dispersed_pop[[i]]) == rownames(adults)[j]), ]
      
    }
  }
  return(dispersed_pop)
}

# Bring them all together for a simulation
simulation <- function(patches,
                       n_per_patch, 
                       coords,
                       n_loci,
                       mate_prob, 
                       bloodmeal_prob, 
                       fecundity, 
                       drive_conversion_prob,
                       daily_survival, 
                       daily_transition,
                       carry_k,
                       sim_days,
                       dispersal_matrix) {
  pop <- ini_pop(patches, n_per_patch, coords, n_loci)
  
  
  patch_sizes <- list()
  stage_distributions <- list()
  
  for (day in 1:sim_days) {
    cat("Day", day, "Completed\n")
    
    
    # Growth with reproduction
    pop <- growth(pop_patches = pop, 
                  mate_prob, 
                  bloodmeal_prob, 
                  fecundity,
                  drive_conversion_prob,
                  n_loci,
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