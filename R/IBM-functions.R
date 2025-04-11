# Functions for running the stephensii IBM

# initialise population
ini_pop <- function(patches, n_per_patch, coords, loci) {
  patches_pop <- list()
  
  for (i in 1:patches) {
    patches_pop[[i]] <- tibble(
      sex = rbinom(n_per_patch[i], 1, 0.5), # Female == 1, random sex
      stage = sample(c("egg", "larva", "pupa", "adult"), n_per_patch[i], replace = TRUE),
      alive = TRUE,
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
  browser()
  updated_pop_patches <- list()
  
  for (i in seq_along(pop_patches)) {
    pop <- pop_patches[[i]]  
    
    male <- pop |> filter(sex == 0, stage == "adult")       # All females
    fem <- pop |> filter(sex == 1, stage == "adult")            # All females
    mated_fem <- rbinom(1, nrow(fem), mate_prob)         # All mated females
    bloodfed_fem <- rbinom(1, mated_fem, bloodmeal_prob)   # Probability that a female finds a blood meal
    
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
    #   allele1 = matrix(sample(c(fem$allele1, fem$allele2), n_offspring * n_loci, replace = TRUE), ncol = n_loci),
    #   allele2 = matrix(sample(c(male$allele1, male$allele2), n_offspring * n_loci, replace = TRUE), ncol = n_loci)
    # )
    
    
    # Genetic inheritance and drive conversion
    # if (n_offspring > 0) { 
    if (n_offspring > 0 && nrow(male$allele1) > 0 && nrow(male$allele2) > 0) {
      dams <- slice_sample(fem, n = bloodfed_fem, replace = FALSE) # sample breeding females
      dams <- slice_sample(dams, n = n_offspring * n_loci, replace = TRUE) |> select(contains("allele")) # assign them randomly to offspring
      sires <- slice_sample(male, n = n_offspring * n_loci, replace = TRUE) |> select(contains("allele")) # randomly assign sires to offspring
      
      sample_allele <- function(allele_tibble){
        select_allele <- rbinom(length(allele_tibble), 1, 0.5) # which dam allele does junior get?
        allele_tibble$allele1*select_allele + (allele_tibble$allele2)*(1-select_allele)
      }
      
      
      dam_allele <- sample_allele(dams)
      sire_allele <- sample_allele(sires)
      
      offspring <- tibble(
        sex = rbinom(n_offspring, 1, 0.5),
        stage = rep("egg", n_offspring),
        alive = TRUE,
        allele1 = dam_allele,
        allele2 = sire_allele
      )
      
      # conversion 
      # converted, or not ### UP TO HERE ###
      converted <- rbinom(length(offspring$allele1), 1, drive_conversion_prob)
      
      
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
dispersal <- function(pop, dispersal_matrix, check = FALSE) {
  #browser()
  patch_indices <- dispersed_pop <- vector(mode = "list", length = nrow(dispersal_matrix))
  
  # get new patch indices for each adult
  for (i in 1:length(pop)) {
    patch <- pop[[i]]
    
    # adults in the patch capable of dispersing
    ad_subset <- patch$stage == "adult" & patch$alive == TRUE
    n_adults <- sum(ad_subset)
    
    # If there are no adults, skip dispersal for this patch
    if (n_adults == 0) next
    
    # Get the dispersal probabilities for this individual according to the dispersal matrix
    dispersal_probs <- dispersal_matrix[i,]
    ad_within_pop_indices <- which(ad_subset)
    new_pop_indices <- sample(1:length(dispersal_probs), size = n_adults, replace = TRUE, prob = dispersal_probs)
    patch_indices[[i]] <- tibble(ad_within_pop_indices, new_pop_indices)
    dispersed_pop[[i]] <- filter(patch, !ad_subset) # move non-adults to same patch in the future
  }
  
  # move individuals to new patches
  for (i in 1:length(pop)) {
    patch <- pop[[i]]
    adults <- patch[patch_indices[[i]]$ad_within_pop_indices, ]
    for (jj in 1:length(pop)){
      ads_jj <- filter(adults, patch_indices[[i]]$new_pop_indices == jj)
      dispersed_pop[[jj]] <- bind_rows(dispersed_pop[[jj]], ads_jj)
    }
  }
  if (check){
    n_pop <- sum(sapply(pop, nrow))
    n_disp <- sum(sapply(dispersed_pop, nrow))
    cat(n_pop, " ", n_disp, "\n")
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

make_dispersal_matrix <- function(coords, lambda, dispersal_frac) {
  # dispersal matrix 
  dist_matrix <- as.matrix(dist(coords, method = "euclidean"))
  
  dispersal_kernel <- exp(-lambda * dist_matrix)
  
  # set the diagonal elements to 0 to prevent self-dispersal
  diag(dispersal_kernel) <- 0
  
  
  # make these rows sum to 1 to get probability of moving to other patch
  # *if* they left. This dispersal matrix gives the probability of the vector
  # vector moving between patches
  rel_dispersal_matrix <- sweep(dispersal_kernel, 1,
                                rowSums(dispersal_kernel), FUN = "/")
  
  # normalise these to have the overall probability of dispersing to that patch,
  # and add back the probability of remaining
  dispersal_matrix <- dispersal_frac * rel_dispersal_matrix +
    (1 - dispersal_frac) * diag(nrow(dispersal_kernel))
  
  return(dispersal_matrix)
}