

# core functions to run Anopheles stephensi population dynamics --------



# Initial population setup ####

# initialise population

ini_pop <- function(patches, n_per_patch, coords, loci, init_frequency) {
  patches_pop <- list()
  
  for (i in 1:patches) {
    patches_pop[[i]] <- tibble(
      sex = rbinom(n_per_patch[i], 1, 0.5), # Female == 1, random sex
      stage = sample(c("egg", "larva", "pupa", "adult"), n_per_patch[i], replace = TRUE),
      allele1 = matrix(rbinom(n = n_per_patch[i] * n_loci, size = 1, prob = init_frequency), ncol = n_loci), # 0 = wild-type, 1 = drive allele
      allele2 = matrix(rbinom(n = n_per_patch[i] * n_loci, size = 1, prob = init_frequency), ncol = n_loci),
      male_allele1 = matrix(NA, nrow = n_per_patch[i], ncol = n_loci),
      male_allele2 = matrix(NA, nrow = n_per_patch[i], ncol = n_loci),
      gdd_accumulated = 0,
      next_oviposition = 0,
      parity1 = 0,
      parity2 = 0,
      parity3 = 0,
      mated = 0,
      fed = 0,
      gravid = 0,
      birth = NA_integer_,
      first_ovip_day = NA_integer_,
      alive = TRUE
    )
    if (length(n_per_patch) != patches) warning("Initial patch population does not equal specified number of patches")
  }
  
  return(patches_pop)
}




# Growth, reproduction and genetic/drive inheritance ####

growth <- function(pop_patches, 
                   bloodmeal_prob, 
                   n_loci,
                   beta,
                   decay,
                   lethal_effect,
                   sterile,
                   sim_days,
                   t_max,
                   t_min,
                   humidty,
                   surface_area,
                   loci_cov_matrix,
                   ldt,
                   mu,
                   sigma_dd) {
    # if (sim_days == 12) browser()
    updated_pop_patches <- list()
    
    for (i in seq_along(pop_patches)) {
      pop <- pop_patches[[i]]  
      
      male <- pop[pop$sex == 0 & pop$stage == "adult", ]  # All males
      fem <- pop[pop$sex == 1 & pop$stage == "adult", ]   # All females
      n.fem <- nrow(fem)
      n.male <- nrow(male)

# Mating 
    if (n.fem > 0 && n.male > 0){
      unmated <- which(fem$mated == 0)
      n.unmated <- length(unmated)
      
      if (n.unmated > 0) {
        realised_mated <- rbinom(n.unmated, 1, prob = (n.male / (beta + n.male)))
        mate_now <- unmated[realised_mated == 1]
        n.mate_now <- length(mate_now)
     
         if (n.mate_now > 0) {
        
        fem$mated[mate_now] <- 1
        selected_male_idx <- sample(n.male, n.mate_now, replace = TRUE)
        selected_male <- male[selected_male_idx,]
        fem$male_allele1[mate_now,] <- selected_male$allele1
        fem$male_allele2[mate_now,] <- selected_male$allele2
        }
      }
      
      
# Blood feeding
      non_fed <- which(fem$mated == 1 & fem$fed == 0)
      n.non_fed <- length(non_fed)
      
      if (n.non_fed > 0) {
        realised_bloodmeal <- rbinom(n.non_fed, 1, bloodmeal_prob)
        fem$fed[non_fed[realised_bloodmeal == 1]] <- 1
      }
      
      #### Update gravid status
      fem$gravid <- ifelse(fem$mated == 1 & fem$fed == 1, 1, 0)
      fem$next_oviposition[fem$gravid == 1] <- fem$next_oviposition[fem$gravid == 1] + 1
      
     
       # estimate clutch sizes and oviposition timing using daily average temperature  
      max_temp <- t_max[i]
      min_temp <- t_min[i]
      daily_temp <- (max_temp+min_temp)/2
      
      delay <- sim_delays(n.fem, daily_temp)
      batch_sizes <- sim_batch_sizes(n.fem)
      
      #### Oviposition conditions 
      cond1 <- as.numeric(fem$next_oviposition >= delay & fem$parity1 == 0 & fem$gravid == 1)
      cond2 <- as.numeric(fem$next_oviposition >= delay & fem$parity1 == 1 & fem$parity2 == 0 & fem$gravid == 1)
      cond3 <- as.numeric(fem$next_oviposition >= delay & fem$parity2 == 1 & fem$parity3 == 0 & fem$gravid == 1)

      # homozygous loci for each female
      homo_loci <- rowSums((fem$allele1 + fem$allele2) == 2)        
      
      # oviposition
      if (sterile) {
        homozygous <- (homo_loci > 0)
        sterile <- as.numeric(!homozygous)
        exp_offspring1 <- cond1 * fem$gravid * batch_sizes * sterile
        exp_offspring2 <- cond2 * fem$gravid * batch_sizes * sterile
        exp_offspring3 <- cond3 * fem$gravid * batch_sizes * sterile
        exp_offspring <- exp_offspring1 + exp_offspring2 + exp_offspring3
      } else {
        exp_offspring1 <- cond1 * fem$gravid * batch_sizes
        exp_offspring2 <- cond2 * fem$gravid * batch_sizes
        exp_offspring3 <- cond3 * fem$gravid * batch_sizes
        exp_offspring <- exp_offspring1 + exp_offspring2 + exp_offspring3
      }
      
      
      fem$parity1[cond1 == 1] <- 1
      fem$first_ovip_day[cond1 == 1 & is.na(fem$first_ovip_day)] <- sim_days
      fem$parity2[cond2 == 1] <- 1
      fem$parity3[cond3 == 1] <- 1
      oviposited <- which((cond1+cond2+cond3) > 0)
      fem$next_oviposition[oviposited] <- 0
      fem$fed[oviposited] <- 0
      fem$gravid[oviposited] <- 0
      
    }   else {
      # If not, set clutch size to 0
      exp_offspring <- rep(0, n.fem)
    }
    
    
    # Offspring generation: Draw the actual number of offspring from a Poisson distribution
    n_offspring <- rpois(n.fem, exp_offspring)
    total_offspring <- sum(n_offspring)
      
    
      if (total_offspring > 0){  
      # Replicate the parents features `n_offspring` times for each offspring, collect only genetic information

        fem_germline <- fem[rep(1:n.fem, n_offspring), c("allele1", "allele2")]
        male_germline <- fem[rep(1:n.fem, n_offspring), c("male_allele1", "male_allele2")]
        
        # Genetic inheritance
        num_loci <- ncol(fem_germline$allele1)
        stopifnot(num_loci == n_loci)
        
        # # random selection of allele, with linkage 
        
        which_allele_fn <- function(n_offspring, num_loci, loci_cov_matrix){
          epsilon <- MASS::mvrnorm(n_offspring, rep(0, num_loci), Sigma = loci_cov_matrix)
          selection_prob <- plogis(epsilon)
          matrix(rbinom(n_offspring * num_loci, 1, selection_prob) == 1,
                 nrow = n_offspring,
                 ncol = num_loci)
        }
        
        # alternative  fucbtion for computational speed
        # 
        # which_allele_fn <- function(n_ind, n_loci, loci_cov_matrix) {
        #   epsilon <- MASS::mvrnorm(n = n_ind,
        #                            mu = rep(0, n_loci),
        #                            Sigma = loci_cov_matrix)
        #   
        #   # alternatively, pass in 'L_loci_cov_matrix', which is computed earlier as:
        #   #   L_loci_cov_matrix <- chol(loci_cov_matrix)
        #   # then inside this function do:
        #   #   z <- matrix(rnorm(n_ind * n_loci), n_ind, n_loci)
        #   #   epsilon <- z %*% L
        #   
        #   selection_prob <- 1 / (1 + exp(-epsilon))
        #   u <- matrix(runif(n_ind * n_loci),
        #               n_ind, n_loci)
        #   u < selection_prob
        # }
    
        which_allele_female <- which_allele_fn(total_offspring, num_loci, loci_cov_matrix) # female gametes
        which_allele_male <- which_allele_fn(total_offspring, num_loci, loci_cov_matrix) # male gametes
      
        #  Determination of offspring features
        offspring <- tibble(
          sex = rbinom(total_offspring, 1, 0.5),
          stage = "egg",
          allele1 = ifelse(which_allele_female,
                         fem_germline$allele1,
                         fem_germline$allele2),
          allele2 = ifelse(which_allele_male,
                         male_germline$male_allele1,
                         male_germline$male_allele2),
          gdd_accumulated = 0,
          next_oviposition = 0,
          parity1 = 0,
          parity2 = 0,
          parity3 = 0,
          mated = 0,
          fed = 0,
          gravid = 0,
          birth = sim_days,
          first_ovip_day = NA_integer_,
          alive = TRUE
        )
      
        # Update pop with offspring & fem population
        pop <- pop[!(pop$sex == 1 & pop$stage == "adult"), ]
        pop <- bind_rows(pop, offspring, fem)
      } else {
        # Update pop with females only
        pop <- pop[!(pop$sex == 1 & pop$stage == "adult"), ]
        pop <- bind_rows(pop, fem)
      }
      
    
      # Genetic load: lethal effect
      
      if (lethal_effect){
        homozygous_lethal <- (pop$allele1 == 1) & (pop$allele2 == 1)
        any_homozygous <- rowSums(homozygous_lethal) > 0
        pop <- filter(pop, !any_homozygous)
        #pop <- pop[pop[!any_homozygous], ]
      }
      
      
    # # Gene Drive architecture (conversion mechanism and inheritance)
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
  
      
    
    # daily humidity  
    daily_humidity <- humidity[i]
      
    
    # Density-dependent survival for aqauatic stages
    count <- sum((pop$stage == "egg") + (pop$stage == "larva") + (pop$stage == "pupa"))
    aquatic_stage_density <- count/surface_area
    
    
    egg_gdd_accumulated <- cal_dd (max_temp, min_temp, ldt["egg"])
    larva_gdd_accumulated <- cal_dd (max_temp, min_temp, ldt["larva"])
    pupa_gdd_accumulated <- cal_dd (max_temp, min_temp, ldt["pupa"])

    
    pop <- pop |>
     mutate(
       gdd_accumulated = case_when(
         stage == "egg"   ~ gdd_accumulated + egg_gdd_accumulated,
         stage == "larva" ~ gdd_accumulated + larva_gdd_accumulated,
         stage == "pupa"  ~ gdd_accumulated + pupa_gdd_accumulated,
         stage == "adult"  ~ 0
       )
     ) 
    
    
   pop <- pop |>
     mutate(
       transition_egg_larva  = stage == "egg"  & rbinom(n(), 1, prob_trans(gdd_accumulated, mu["egg"], sigma_dd["egg"])) == 1,
       transition_larva_pupa = stage == "larva" & rbinom(n(), 1, prob_trans(gdd_accumulated, mu["larva"], sigma_dd["larva"])) == 1,
       transition_pupa_adult = stage == "pupa"  & rbinom(n(), 1, prob_trans(gdd_accumulated, mu["pupa"], sigma_dd["pupa"])) == 1,
     ) |>
     mutate(
       stage = case_when(
         transition_egg_larva  ~ "larva",
         transition_larva_pupa ~ "pupa",
         transition_pupa_adult ~ "adult",
         TRUE ~ stage
       ),
       gdd_accumulated = case_when(
         transition_egg_larva  ~ 0,
         transition_larva_pupa ~ 0,
         transition_pupa_adult ~ 0,
         TRUE ~ gdd_accumulated
       )
     ) |>
     select(-starts_with("transition_"))
   
   
   # Density dependent survival adjusted to temperature (for aquatic stage) and humidty (for adult stage)
  
     pop <- pop |> mutate(
       alive = case_when(
         stage == "egg" ~ rbinom(n(), 1, das_temp_dens_As(daily_temp, aquatic_stage_density)),
         stage == "larva" ~ rbinom(n(), 1, das_temp_dens_As(daily_temp, aquatic_stage_density)),
         stage == "pupa" ~ rbinom(n(), 1, das_temp_dens_As(daily_temp, aquatic_stage_density)),
         stage == "adult" ~ rbinom(n(), 1, ds_temp_humid_As(daily_temp, daily_humidity, species = "An. stephensi")),
         TRUE ~ NA_integer_
       ),
       alive = alive == 1
   )
   

    pop <- pop[pop$alive,]
    
    updated_pop_patches[[i]] <- pop
  }
  return(updated_pop_patches)
 }


# Dispersal ####

#### Metapopulation dispersal function ####
  
  meta_dispersal <- function(pop, dispersal_matrix, check = FALSE) {
    
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
      dispersed_pop[[i]] <- patch[!ad_subset, ]
    }
    
    # move individuals to new patches
    for (i in 1:length(pop)) {
      patch <- pop[[i]]
      adults <- patch[patch_indices[[i]]$ad_within_pop_indices, ]
      for (jj in 1:length(pop)){
        ads_jj <- adults[patch_indices[[i]]$new_pop_indices == jj, ]
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



#### Stepping stone dispersal (one-dimensional discrete space) ####
  

    ss_dispersal <- function(pop_patches, dispersal_prop) {
      n_patches <- length(pop_patches)
      dispersed_pop <- vector("list", n_patches)
      
      # Initialize empty population for each patch
      for (i in seq_len(n_patches)) {
        dispersed_pop[[i]] <- pop_patches[[i]][0, ]
      }
      
      for (i in seq_len(n_patches)) {
        current_patch <- pop_patches[[i]]
        
        if (nrow(current_patch) == 0) next
        
        # Split adults and non-adults
        adults <- current_patch[current_patch$stage == "adult", ]
        # can_disperse <- rbinom(nrow(current_patch[current_patch$stage == "adult", ]), 1, dispersal_prop)
        non_adults <- current_patch[current_patch$stage != "adult", ]
        
        # Determine which adults disperse
        ready_to_disperse <- rbinom(nrow(adults), 1, dispersal_prop)
        dispersing_adults <- adults[ready_to_disperse == 1, ]
        # dispersing_adults <- current_patch[can_disperse == 1, ]
        staying_adults <- adults[ready_to_disperse == 0, ]
        
        # Disperse adults to i-1 or i+1
        if (nrow(dispersing_adults) > 0) {
          directions <- sample(c(-1, 1), nrow(dispersing_adults), replace = TRUE)
          target_patch <- i + directions
          target_patch <- pmin(pmax(target_patch, 1), n_patches)  # keep within boundaries i.e. patch 1 and patch "n"
          
          for (j in seq_along(target_patch)) {
            dispersed_pop[[target_patch[j]]] <- bind_rows(
              dispersed_pop[[target_patch[j]]],
              dispersing_adults[j, ]
            )
          }
        }
        
        # Add stayers (non-dispersing adults and non-adults) to current patch
        dispersed_pop[[i]] <- bind_rows(
          dispersed_pop[[i]], 
          staying_adults,
          non_adults
        )
      }
      
      return(dispersed_pop)
    }

  

# function to run simulation ####
run_model <- function(patches,
                       n_per_patch, 
                       coords,
                       n_loci,
                       init_frequency,
                       bloodmeal_prob, 
                       beta,
                       decay,
                       lethal_effect,
                       sterile,
                       sim_days,
                       stepping_stone_model,
                       dispersal_matrix,
                       t_max,
                       t_min,
                       humidty,
                       surface_area,
                       ldt,
                       mu,
                       sigma_dd) {
  
  pop <- ini_pop(patches, n_per_patch, coords, n_loci, init_frequency)
  
  patch_sizes <- list()
  allele_frequency <- list()
  #spread_rate <- list()
  #generation_time list()
  
  
  for (day in 1:sim_days) {
    #if (day == 45) browser()
    cat("Day", day, "Underway \n")
    # Growth with reproduction
    pop <- growth(pop_patches = pop,
                  bloodmeal_prob, 
                  n_loci,
                  beta,
                  decay,
                  lethal_effect,
                  sterile,
                  sim_days = day,
                  t_max,
                  t_min,
                  humidty,
                  surface_area,
                  loci_cov_matrix = l.cov.mat,
                  ldt,
                  mu,
                  sigma_dd)
    
    # Dispersal
    if(stepping_stone_model) {
      pop <- ss_dispersal(pop, dispersal_prop)
    } else {
      pop <- meta_dispersal(pop, dispersal_matrix, check = FALSE)
    }
    
    # Track daily population sizes per patch
  
    patch_sizes[[day]] <- do.call(rbind, lapply(seq_along(pop), function(patch_id) {
      data.frame(
        day = day,
        patch = patch_id,
        pop_size = nrow(pop[[patch_id]])
      )
    }))
    patch_sizes_df <- do.call(rbind, patch_sizes)
    
    
    # Track daily allele frequency per patch

    allele_frequency[[day]] <- do.call(rbind, lapply(seq_along(pop), function(patch_id) {
      patch_pop <- pop[[patch_id]]
      if (nrow(patch_pop) > 0) {
        deleterious <- sum(patch_pop$allele1 == 1) + sum(patch_pop$allele2 == 1)
        total <- 2 * nrow(patch_pop) * ncol(patch_pop$allele1)
        wild_type <- total - deleterious
        freq <- deleterious / total
      } else {
        deleterious <- 0
        total <- 0
        wild_type <- 0
        freq <- 0
      }
      data.frame(
        patch = patch_id,
        wild = wild_type,
        lethal = deleterious,
        total = total,
        freq = freq,
        day = day
      )
    }
    ))
    allele_frequency_df <- do.call(rbind, allele_frequency)
  }
  
  # track spread or invasion rate
  
  # Return the collected data
  list(
    patch_sizes = patch_sizes_df,
    allele_frequency = allele_frequency_df,
    # spread_rate <-
    final_pop = pop
  )
}


