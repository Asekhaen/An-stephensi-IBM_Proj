

# core functions to run Anopheles stephensi population dynamics --------

# Initial population setup ####


ini_pop <- function(patches, n_per_patch, coords, n_loci) { 
  patches_pop <- list()
  
  for (i in 1:patches) {
    patches_pop[[i]] <- tibble(
      stage = sample(stages, n_per_patch[i], replace = TRUE),
      # chromosome1 = matrix(rbinom(n = n_per_patch[i] * n_loci, size = 0, prob = init_frequency), ncol = n_loci), # 0 = wild-type, 1 = drive allele
      # chromosome2 = matrix(rbinom(n = n_per_patch[i] * n_loci, size = 0, prob = init_frequency), ncol = n_loci),
      allo1 = make_allosome("X", n_per_patch[i]),
      allo2 = make_allosome(c("X", "Y"), n_per_patch[i]),
      sex = if_else(
        allo1[,1] == "X" & allo2[,1] == "X",
        "female",
        "male"),
      autosome1 = matrix(0, nrow = n_per_patch[i], ncol = n_loci),    # 0 = wild type, 1 = drive, 2 = resistance
      autosome2 = matrix(0, nrow =n_per_patch[i], ncol = n_loci),
      male_allo1 = matrix(NA_character_, nrow = n_per_patch[i], ncol = 2),
      male_allo2 = matrix(NA_character_, nrow = n_per_patch[i], ncol = 2),
      male_autosome1 = matrix(NA, nrow = n_per_patch[i], ncol = n_loci),
      male_autosome2 = matrix(NA, nrow = n_per_patch[i], ncol = n_loci),
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
     # if (sim_days == 3) browser()
     #browser()
    updated_pop_patches <- list()
    
    for (i in seq_along(pop_patches)) {
      pop <- pop_patches[[i]]  
      
      male <- pop[pop$sex == "male" & pop$stage == "adult", ]  # All males
      fem <- pop[pop$sex == "female" & pop$stage == "adult", ]   # All females
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
        fem$male_autosome1[mate_now,] <- selected_male$autosome1
        fem$male_autosome2[mate_now,] <- selected_male$autosome2
        fem$male_allo1[mate_now,] <- selected_male$allo1
        fem$male_allo2[mate_now,] <- selected_male$allo2
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
      
     
      # estimate egg clutch size and timing of oviposition using daily average temperature  
      max_temp <- temp_max[i]
      min_temp <- temp_min[i]
      daily_temp <- (max_temp+min_temp)/2
      
      delay <- sim_delays(n.fem, daily_temp)
      batch_sizes <- sim_batch_sizes(n.fem)
      
      #### Oviposition conditions. This simulates the intervals between when the 
      #### female becomes gravid and oviposition
      cond1 <- as.numeric(fem$next_oviposition >= delay & fem$parity1 == 0 & fem$gravid == 1)
      cond2 <- as.numeric(fem$next_oviposition >= delay & fem$parity1 == 1 & fem$parity2 == 0 & fem$gravid == 1)
      cond3 <- as.numeric(fem$next_oviposition >= delay & fem$parity2 == 1 & fem$parity3 == 0 & fem$gravid == 1)


      
      # homo_loci <- rowSums(
      #   (fem$autosome1 + fem$autosome2) == 2)
      
      
      wt_homozygous  <- fem$autosome1 == 0 & fem$autosome2 == 0
      drive_homozygous <- fem$autosome1 == 1 & fem$autosome2 == 1  # sterile
      res_homozygous <- fem$autosome1 == 2 & fem$autosome2 == 2   # sterile
      drive_wt <- fem$autosome1 == 0 & fem$autosome2 == 1 | fem$autosome1 == 1 & fem$autosome2 == 0
      drive_res <- fem$autosome1 == 1 & fem$autosome2 == 2 | fem$autosome1 == 2 & fem$autosome2 == 1   # sterile
      res_wt <- fem$autosome1 == 0 & fem$autosome2 == 2 | fem$autosome1 == 2 & fem$autosome2 == 0
      
      disrupted_loci <- drive_homozygous + res_homozygous + drive_res


      # oviposition (with effect of deleterious allele on fitness: sterility)
      
      if (sterile) {
        sterile_loci <- (disrupted_loci > 0)
        sterility <- as.integer(rowSums(sterile_loci) != ncol(sterile_loci))
        exp_offspring1 <- cond1 * fem$gravid * batch_sizes * sterility
        exp_offspring2 <- cond2 * fem$gravid * batch_sizes * sterility
        exp_offspring3 <- cond3 * fem$gravid * batch_sizes * sterility
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
      # If not, set clutch size to 0. i.e no offspring
      exp_offspring <- rep(0, n.fem)
    }
    
    # I changed NA to "zeros" because the operation produced NAs from none mated 
    # individuals that are still part of the female population 
    exp_offspring <- replace(exp_offspring, is.na(exp_offspring), 0) 
    
    # Offspring generation: Draw the actual number of offspring from a Poisson distribution
    n_offspring <- rpois(n.fem, exp_offspring)
    total_offspring <- sum(n_offspring, na.rm = TRUE)
      
    
      if (total_offspring > 0){  
      # Replicate the parents features `n_offspring` times for each offspring, collect only genetic information

        fem_germline <- fem[rep(1:n.fem, n_offspring), c("autosome1", "autosome2")]
        male_germline <- fem[rep(1:n.fem, n_offspring), c("male_autosome1", "male_autosome2")]
        
        fem_allo <- fem[rep(1:n.fem, n_offspring), c("allo1", "allo2")]
        male_allo <- fem[rep(1:n.fem, n_offspring), c("male_allo1", "male_allo2")]
        
        
        
        # Genetic inheritance
        num_loci <- ncol(fem_germline$autosome1)
        stopifnot(num_loci == n_loci)
        
        # # random selection for linked loci 
        
        # which_allele_fn <- function(n_offspring, num_loci, loci_cov_matrix){
        #   epsilon <- MASS::mvrnorm(n_offspring, rep(0, num_loci), Sigma = loci_cov_matrix)
        #   selection_prob <- plogis(epsilon)
        #   matrix(rbinom(n_offspring * num_loci, 1, selection_prob) == 1,
        #          nrow = n_offspring,
        #          ncol = num_loci)
        # }
        
        # alternative  function for computational speed
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
    
        # which_allele_female <- which_allele_fn(total_offspring, num_loci, loci_cov_matrix) # female gametes
        # which_allele_male <- which_allele_fn(total_offspring, num_loci, loci_cov_matrix) # male gametes
      
        
        
        which_autosome  <- matrix(rep(rbinom(total_offspring, 1, 0.5), n_loci),
                                    nrow = total_offspring, ncol = n_loci)
        which_autosome_mate <- matrix(rep(rbinom(total_offspring, 1, 0.5), n_loci),
                                    nrow = total_offspring, ncol = n_loci)
        
        
        
        which_allosome  <- matrix(rep(rbinom(total_offspring, 1, 0.5), 2),
                                  nrow = total_offspring, ncol = 2)
        which_allosome_mate <- matrix(rep(rbinom(total_offspring, 1, 0.5), 2),
                                      nrow = total_offspring, ncol = 2)
        
        
        #  Determination of offspring features
        offspring <- tibble(
          stage = "egg",
          allo1 = ifelse(which_allosome,
                         fem_allo$allo1,
                         fem_allo$allo2),           #if (rbinom(1, 1, 0.5) == 1) fem_allo$allo1 else fem_allo$allo2,
          allo2 = ifelse(which_allosome_mate,
                         male_allo$male_allo1,
                         male_allo$male_allo2),     #if (rbinom(1, 1, 0.5) == 1) male_allo$male_allo1 else male_allo$male_allo2,
          sex = ifelse(allo1[,1] == "X" & allo2[,1] == "X", "female", "male"),
          autosome1 = ifelse(which_autosome,
                             fem_germline$autosome1,
                             fem_germline$autosome2), #if (rbinom(1, 1, 0.5) == 1) fem_germline$autosome1 else fem_germline$autosome2,
          autosome2 =  ifelse(which_autosome_mate,
                              male_germline$male_autosome1,
                              male_germline$male_autosome2),#if (rbinom(1, 1, 0.5) == 1) male_germline$male_autosome1 else male_germline$male_autosome2,
          male_autosome1 = matrix(NA, ncol = n_loci),
          male_autosome2 = matrix(NA, ncol = n_loci),
          male_allo1 = matrix(NA_character_, ncol = 2),
          male_allo2 = matrix(NA_character_, ncol = 2),
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
        pop <- pop[!(pop$sex == "female" & pop$stage == "adult"), ]
        pop <- bind_rows(pop, offspring, fem)
      } else {
        # Update pop with females only if no oviposition happened 
        pop <- pop[!(pop$sex == "female" & pop$stage == "adult"), ]
        pop <- bind_rows(pop, fem)
      }
      
      # # effect of deleterious allele on fitness: lethal effect
      # 
      # if (lethal_effect){
      #   homozygous_lethal <- (pop$autosome1 == 1) & (pop$autosome2 == 1)
      #   any_homozygous <- rowSums(homozygous_lethal) > 0
      #   pop <- filter(pop, !any_homozygous)
      #   #pop <- pop[pop[!any_homozygous], ]
      # }
      # 
      
    # # Gene Drive architecture (conversion mechanism and inheritance)
    # 
    #     drive_conversion <- function(parent, prob1, prob2) {
    #       if (any(is.na(parent$chromosome1)) || any(is.na(parent$chromosome2))) {
    #           warning("NA detected in allele input!")
    #         }
    #       heterozygous <- (parent$chromosome1 + parent$chromosome2) == 1
    #       # Drive conversion 95% conversion rate
    #       converted <- rbinom(length(parent$chromosome1), 1, prob1) # drive conversion at each locus
    #       conv_event <- converted*heterozygous # conversion event?
    #       parent$chromosome1[parent$chromosome1 == 0 & conv_event == 1] <- 1 # successful conversions
    #       parent$chromosome2[parent$chromosome2 == 0 & conv_event == 1] <- 1
    # 
    #       # Resistance development (0 → 2)
    #       failed_conv <- heterozygous & conv_event == 0
    #       resistance_event <- rbinom(length(parent$chromosome1), 1, prob2)
    #       parent$chromosome1[parent$chromosome1 == 0 & failed_conv & resistance_event == 1] <- 2  # Thoughts/To do: individuals that did not develop resistance, yet heterozygous can be can be designated as those with functional resistance and resistant to future Cas9 cutting
    #       parent$chromosome2[parent$chromosome2 == 0 & failed_conv & resistance_event == 1] <- 2
    # 
    #       return(parent)
    #     }
    #     
    #     
    #     fem_germline <- drive_conversion(fem_germline, conversion_prob, resistance_prob)   # dams with drive converted germ line
    #     male_germline <- drive_conversion(male_germline, conversion_prob, resistance_prob) # sires with drive converted germ line
  
      
  
    # stage development using growth-degree day accumulation
 
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
   
   
   # Density-dependent survival adjusted to temperature + density (aquatic stage) and temperature + humidity (adult stage)
   
   # density of the auqtic stages (eggs, larvae, and pupae) 
   count <- sum((pop$stage == "egg") + (pop$stage == "larva") + (pop$stage == "pupa"))
   aquatic_stage_density <- count/surface_area
   # daily humidity  
   daily_humidity <- humidity[i]
  
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


  dispersal <- function(pop, dispersal_type, check = FALSE) {
    
    patch_indices <- dispersed_pop <- vector(mode = "list", length = nrow(dispersal_type))
    
    # get new patch indices for each adult
    for (i in 1:length(pop)) { #define this first off)
      patch <- pop[[i]]
      
      # adults in the patch capable of dispersing
      dispersal_ready <- patch$stage == "adult" & patch$alive == TRUE
      n_dispersal_ready <- sum(dispersal_ready)
      
      # If there are no adults, skip dispersal for this patch
      if (n_dispersal_ready == 0) next
      
      # Get the dispersal probabilities for this individual according to the dispersal matrix
      dispersal_probs <- dispersal_type[i,]
      dispersers <- which(dispersal_ready)
      new_pop_indices <- sample(1:length(dispersal_probs), size = n_dispersal_ready, replace = TRUE, prob = dispersal_probs)
      patch_indices[[i]] <- tibble(dispersers, new_pop_indices)
      dispersed_pop[[i]] <- patch[!dispersal_ready, ] 
    }
    
    # move individuals to new patches
    for (i in 1:length(pop)) {
      patch <- pop[[i]]
      adults <- patch[patch_indices[[i]]$dispersers, ]
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
                       dispersal_type,
                       t_max,
                       t_min,
                       humidty,
                       surface_area,
                       ldt,
                       mu,
                       sigma_dd) {
  
  pop <- ini_pop(patches, n_per_patch, coords, n_loci)
  
  patch_sizes <- list()
  #colonisation_rate_list <- list()
  allele_frequency <- list()
  # invasion_speed <- list()
  
  
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
  
      pop <- dispersal(pop, dispersal_type, check = FALSE)
    
    # Track daily population sizes, patch occupancy rates, etc.
    patch_sizes[[day]] <- tibble(
      day = day,
      patch = seq_along(pop),
      pop_size = sapply(pop, nrow),
      patch_occupied = sum(pop_size > 0),
      unoccupied = length(patch) - patch_occupied,
      occupancy_rate = patch_occupied/length(patch)
    )
    patch_sizes_df <- bind_rows(patch_sizes)

    # # Track overall allele frequency and allele frequency per locus
    # allele_frequency[[day]] <-  lapply(seq_along(pop), function(patch_id) {
    #   patch_pop <- pop[[patch_id]]
    #   loci_n  <- ncol(patch_pop$autosome1)  
    #   n_ind   <- nrow(patch_pop$autosome1)  
    #   total_allele_overall <- 2 * n_ind * loci_n
    # 
    #   overall <- tibble(
    #     day        = day,
    #     patch      = patch_id,
    #     total      = total_allele_overall,
    #     deleterious= sum(patch_pop$autosome1 == 1) + sum(patch_pop$autosome2 == 1),
    #     wild       = total_allele_overall - deleterious,
    #     freq       = ifelse(total_allele_overall == 0, 0, deleterious / total_allele_overall)
    #   )
    # })
    #   allele_frequency_df <- bind_rows(allele_frequency)
    
    
    
    # # Track allele frequency for overall allele and per locus.... in progress
    # allele_frequency[[day]] <-  lapply(seq_along(pop), function(patch_id) {
    #   patch_pop <- pop[[patch_id]]
    #   loci_n  <- ncol(patch_pop$autosome1)  
    #   n_ind   <- nrow(patch_pop$autosome1)  
    #   total_allele_overall <- 2 * n_ind * loci_n
    #   total_allele_locus <- 2 * n_ind
    #   
    #   overall <- tibble(
    #     day        = day,
    #     patch      = patch_id,
    #     total      = total_allele_overall,
    #     deleterious= sum(patch_pop$autosome1 == 1) + sum(patch_pop$autosome2 == 1),
    #     wild       = total_allele_overall - deleterious,
    #     freq       = deleterious / total_allele_overall
    #   )
    #   
    #   per_locus <- tibble(
    #     day        = day,
    #     patch      = patch_id,
    #     total      = total_allele_locus,
    #     loci = seq_len(loci_n),
    #     deleterious= colSums(patch_pop$autosome1 == 1) + colSums(patch_pop$autosome2 == 1),
    #     wild       = total_allele_locus - deleterious,
    #     freq       = deleterious / total_allele_locus
    #   )
    #     list(per_locus_freq = per_locus, overall_freq = overall)
    # })
    # 
    # overall_df   <- map_dfr(allele_frequency[[day]], "overall_freq")
    # per_locus_df <- map_dfr(allele_frequency[[day]], "per_locus_freq")
   
  
    # invasion speed <- list()
    

  }
  

  # Return the collected data
  list(
    pop_sizes = patch_sizes_df,
    # patch_colonisation_rate = colonisation_df,
    # allele_freq = allele_frequency_df,
    # per_locus_output = per_locus_df,
    # overall_loci = overall_df,
    # invasion_speed <- invasion_speed_df
    final_pop = pop
  )
}


