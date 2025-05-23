# Functions for running the stephensi IBM


# create coordinates for the patches/locations 
coords <- as.data.frame(100 * matrix(runif(patches * 2), ncol = 2))
colnames(coords) <- c("x","y")


#### Initialise population ####

ini_pop <- function(patches, n_per_patch, coords, loci) {
  patches_pop <- list()
  
  for (i in 1:patches) {
    patches_pop[[i]] <- tibble(
      sex = rbinom(n_per_patch[i], 1, 0.5), # Female == 1, random sex
      stage = sample(c("egg", "larva", "pupa", "adult"), n_per_patch[i], replace = TRUE),
      alive = TRUE,
      allele1 = matrix(sample(c(0, 1), n_per_patch[i] * n_loci, replace = TRUE, prob = c(0.95,0.05)), ncol = n_loci,), # 0 = wild-type, 1 = drive allele
      allele2 = matrix(sample(c(0, 1), n_per_patch[i] * n_loci, replace = TRUE, prob = c(1,0)), ncol = n_loci)
    )
    if (length(n_per_patch) != patches) warning("Initial patch population does not equal specified number of patches")
  }
  
  return(patches_pop)
}


#### Growth, reproduction and drive inheritance ####
growth <- function(pop_patches, 
                   bloodmeal_prob, 
                   fecundity,
                   conversion_prob,
                   resistance_prob,
                   n_loci,
                   daily_survival,
                   daily_transition,
                   alpha,
                   beta,
                   decay,
                   effect_size_factor,
                   sim_days,
                   daily_temp,
                   sigma) {
  
  #if (sim_days == 15) browser()
  
  updated_pop_patches <- list()
  
  for (i in seq_along(pop_patches)) {
    pop <- pop_patches[[i]]  

    male <- pop |> filter(sex == 0, stage == "adult")    # All males
    fem <- pop |> filter(sex == 1, stage == "adult")     # All females

    offspring_list <- list()
    if (nrow(fem) > 0 && nrow(male) > 0)
      # Loop through each female to simulate mating, bloodfeeding, and egg/ofspring production
      for (j in 1:nrow(fem)) {
        
        # select a mate
        selected_male_index <- sample(nrow(male), 1, replace = TRUE)
        selected_male <- male[selected_male_index, ]
        
        # Bernoulli trial for mating and feeding (1 if mated, 0 if not)

        if (rbinom(1, 1, (nrow(male)/(beta + nrow(male)))) == 1 && rbinom(1, 1, bloodmeal_prob) == 1) {   #  (nrow(male)/(beta + nrow(male)))) is the mating probability which increases as male population increases (North and Godfray; Malar J (2018) 17:140) 
         
          # If bloodfed, calculate expected offspring for this female
          # exp_offspring <- fecundity
          
          
          # Effect size of genetic load: additive effect from 0 to 1 as the 
          # number loci homozygous for the lethal gene increases
           
          no_homo_loci <- sum(fem[j,]$allele1 == 1 & fem[j,]$allele2 == 1)
          exp_offspring <- round(fecundity * exp(-effect_size_factor * no_homo_loci))
          
          # Draw the actual number of offspring from a Poisson distribution
          n_offspring <- rpois(1, exp_offspring)
          
        } else {
          # If not, set offspring count to 0
          n_offspring <- 0
        }
        
        
        if (n_offspring > 0) {
          dams_data <- fem[j, ]             # Extract the mother's features
          sires_data <- selected_male       # Extract the father's features
          
          # Replicate the parents features `n_offspring` times for each offspring
          fem_germline <- dams_data[rep(1, n_offspring), ]
          male_germline <- sires_data[rep(1, n_offspring), ]
          
          # select the genetic characteristcs 
          
          male_germline <- male_germline |> select(contains("allele"))
          fem_germline <- fem_germline |> select(contains("allele"))
          
    # # Genetic inheritance and Drive conversion
    # 
    #     drive_conversion <- function(parent, prob1, prob2) {
    #       if (any(is.na(parent$allele1)) || any(is.na(parent$allele2))) {
    #           warning("NA detected in allele input!")
    #         }
    #       heterozygous <- (parent$allele1 + parent$allele2) == 1
    #       # Drive conversion 95% conversion rate
    #       converted <- rbinom(length(parent$allele1), 1, prob1) # drive conversion at each locus
    #       conv_event <- converted*heterozygous # conversion event?
    #       parent$allele1[parent$allele1 == 0 & conv_event == 1] <- 1 # successful conversions
    #       parent$allele2[parent$allele2 == 0 & conv_event == 1] <- 1
    # 
    #       # Resistance development (0 → 2)
    #       failed_conv <- heterozygous & conv_event == 0
    #       resistance_event <- rbinom(length(parent$allele1), 1, prob2)
    #       parent$allele1[parent$allele1 == 0 & failed_conv & resistance_event == 1] <- 2  # Thoughts/To do: individuals that did not develop resistance, yet heterozygous can be can be designated as those with functional resistance and resistant to future Cas9 cutting
    #       parent$allele2[parent$allele2 == 0 & failed_conv & resistance_event == 1] <- 2
    # 
    #       return(parent)
    #     }
    #     
    #     
    #     fem_germline <- drive_conversion(fem_germline, conversion_prob, resistance_prob)   # dams with drive converted germ line
    #     male_germline <- drive_conversion(male_germline, conversion_prob, resistance_prob) # sires with drive converted germ line
      
      
        # Genetic inheritance
        total_offspring <- nrow(fem_germline)
        stopifnot(total_offspring == n_offspring)
        
        num_loci <- ncol(fem_germline$allele1)
        stopifnot(num_loci == n_loci)
        
        
        # a function to place loci at random on the genome (of size = 1)
        # also takes exponential decay and variance to produce variance-covariance matrix
        
        place_loci_mat <- function(loci, genome.size = 1, var = 0.5, decay, mean_prob = 0.5){
          loci_positions <- runif(loci, max = genome.size)
          loci_dist_matrix <- as.matrix(dist(loci_positions))^2 
          loci_cov_matrix <- var*exp(-decay*loci_dist_matrix)
          epsilon <- MASS::mvrnorm(1, rep(0, loci), Sigma = loci_cov_matrix)
          selection_prob <- plogis(qlogis(mean_prob) + epsilon)
          
          return(list(loci_positions = loci_positions, 
                      loci_dist_matrix = loci_dist_matrix, 
                      loci_cov_matrix = loci_cov_matrix,
                      epsilon = epsilon,
                      selection_prob = selection_prob))
        }
        
        
        linkage_selection_prob <- place_loci_mat(loci = n_loci, 
                                                 decay = decay)
        
        # random selection of allele
        
        which_allele <- matrix(rbinom(total_offspring * num_loci, 1, linkage_selection_prob$selection_prob) == 0,
                               nrow = total_offspring,
                               ncol = num_loci)
        
        
        #  Determination of offspring features
        offspring <- tibble(
          sex = rbinom(n_offspring, 1, 0.5),
          stage = "egg",
          alive = TRUE,
          allele1 = ifelse(which_allele,
                           fem_germline$allele1,
                           fem_germline$allele2),
          allele2 = ifelse(which_allele,
                           male_germline$allele1,
                           male_germline$allele2)
        )

    
    # # Gene drive lethal effect: mortality increases to 100% as homozygous locus reaches max.
   
      # homozygous_drive <- (offspring$allele1 == 1) & (offspring$allele2 == 1)
      # num_drive_homozygous <- rowSums(homozygous_drive)
      # max_homozygous <- ncol(offspring$allele1)
      # survival_prob <- 1 - (num_drive_homozygous/max_homozygous)
      # offspring$alive <- rbinom(n = length(survival_prob), size = 1, prob = survival_prob) == 1
      # offspring <- filter(offspring, alive)
      
    # Add offspring to main population
 
        offspring_list[[length(offspring_list) + 1]] <- offspring
      }
  }
  
  # Combine all offspring data frames into one final data frame
  offspring_df <- bind_rows(offspring_list)
  
  # add offspring to the main population. 
  pop <- bind_rows(pop, offspring_df)
  
  
    # Temperature-adjusted survival for larval and adult population
    
    temp_effect <- exp(-((daily_temp - 25)^2) / (2 * sigma^2)) # Gaussian process (after Beck-Johnson et al., 2013). 
    temp_adjusted_survival <- daily_survival[c("adult","larva")] * temp_effect
    
    
    # Density-dependent survival for larval stage
    larva_count <- sum(pop$stage == "larva")
    density_dependence <- 1/(1 + (alpha*larva_count))
    density_dependent_survival <- temp_adjusted_survival["larva"] * density_dependence
    

    pop <- pop |> mutate(
        alive = case_when(
          stage == "egg" ~ rbinom(n(), 1, daily_survival["egg"]),
          stage == "larva" ~ rbinom(n(), 1, density_dependent_survival),
          stage == "pupa" ~ rbinom(n(), 1, daily_survival["pupa"]),
          stage == "adult" ~ rbinom(n(), 1, temp_adjusted_survival["adult"]),
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
          TRUE ~ stage
        )
      )
    
    pop <- filter(pop, alive)
    
    updated_pop_patches[[i]] <- pop
  }
  return(updated_pop_patches)
}



#### Make dispersal matrix ####
make_dispersal_matrix <- function(coords, lambda, dispersal_frac) {
  # dispersal matrix 
  dist_matrix <- as.matrix(dist(coords, method = "euclidean"))
  
  #exponential dispersal kernel
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


# create a dispersal matrix
dispersal_matrix <- make_dispersal_matrix(coords = coords, 
                                          lambda = lambda, 
                                          dispersal_frac = dispersal_frac)

#### Dispersal function ####
dispersal <- function(pop, dispersal_matrix, check = FALSE) {
 
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



#### Simulation: Bring them all together ####
simulation <- function(patches,
                       n_per_patch, 
                       coords,
                       n_loci, 
                       bloodmeal_prob, 
                       fecundity, 
                       conversion_prob,
                       resistance_prob,
                       daily_survival, 
                       daily_transition,
                       alpha,
                       beta,
                       decay,
                       effect_size_factor,
                       sim_days,
                       dispersal_matrix,
                       daily_temp,
                       sigma) {
  pop <- ini_pop(patches, n_per_patch, coords, n_loci)
  
  patch_sizes <- list()
  stage_distributions <- list()
  for (day in 1:sim_days) {
    #if (day == 6) browser()
    cat("Day", day, "Completed\n")

    # Growth with reproduction
    pop <- growth(pop_patches = pop,
                  bloodmeal_prob, 
                  fecundity,
                  conversion_prob,
                  resistance_prob,
                  n_loci,
                  daily_survival,
                  daily_transition,
                  alpha,
                  beta,
                  decay,
                  effect_size_factor,
                  sim_days = day,
                  daily_temp = temp[day],
                  sigma)
    
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


