

# core functions to run Anopheles stephensi population dynamics --------

# Initial population setup ####

ini_pop <- function(patches, 
                    n_per_patch, 
                    coords, 
                    n_loci, 
                    prob_wildtype1,
                    prob_wildtype2) {
  patches_pop <- list()
  for (i in 1:patches) {
    patches_pop[[i]] <- tibble(
      stage = rep("adult", n_per_patch[i]),
      allo1 = make_allosome("X", n_per_patch[i]),
      allo2 = make_allosome(c("X", "Y"), n_per_patch[i]),
      sex = if_else(allo1[,1] == "X" & allo2[,1] == "X", "female", "male"),
      autosome1 = matrix(sample(
        c(0,1),                             # polymorphic locus: 0 = wild-type 1, 1 = wild-type 2
        n_loci * n_per_patch[i],
        replace = TRUE,
        prob = c(prob_wildtype1, prob_wildtype2)), 
        ncol = n_loci),
      autosome2 =  matrix(sample(
        c(0,1),                             # polymorphic locus: 0 = wild-type 1, 1 = wild-type 2
        n_loci *n_per_patch[i],
        replace = TRUE,
        prob = c(prob_wildtype1, prob_wildtype2)), 
        ncol = n_loci),
      male_allo1 = matrix(NA_character_, nrow = n_per_patch[i], ncol = 1),
      male_allo2 = matrix(NA_character_, nrow = n_per_patch[i], ncol = 1),
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
    if (length(n_per_patch) != patches) warning(
      "Initial patch population does not equal specified number of patches")
  }

  return(patches_pop)
}

# this function introduces individuals with drives into the population
introduce_drive <- function(initial_pop, per_release, n_loci) {
  drive_adult <- per_release * initial_pop
  which_chr =  matrix(
    rbinom(drive_adult * n_loci, 1, 0.5),
    nrow = drive_adult,
    ncol = n_loci
  )
  autosome1 = matrix(rbinom(drive_adult, 0, 0.5), n_loci,
                     nrow = drive_adult, ncol = n_loci)
  autosome2 = matrix(rbinom(drive_adult, 0, 0.5), n_loci,
                     nrow = drive_adult, ncol = n_loci)
  autosome1 <- which_chr * 2        # 2 denotes drive allele
  autosome2 <- (1 - which_chr) * 2  # I did this to ensure heterozygous on either chromosome only  

  tibble(
    stage = rep("adult", drive_adult),
    allo1 = make_allosome("X", drive_adult),
    allo2 = make_allosome(c("X", "Y"), drive_adult),
    sex = if_else(allo1[,1] == "X" & allo2[,1] == "X", "female", "male"),
    autosome1 = autosome1,
    autosome2 = autosome2,
    male_allo1 = matrix(NA_character_, nrow = drive_adult),
    male_allo2 = matrix(NA_character_, nrow = drive_adult),
    male_autosome1 = matrix(NA, nrow = drive_adult, ncol = n_loci),
    male_autosome2 = matrix(NA, nrow = drive_adult, ncol = n_loci),
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
}


# Growth, reproduction and genetic/drive inheritance ####

growth <- function(pop_patches, 
                   bloodmeal_prob, 
                   n_loci,
                   beta,
                   decay,
                   recomb,
                   # lethal_effect,
                   drive_effect,
                   drive_type,
                   prob1,
                   prob2,
                   homing_prob,
                   shred_prob,
                   sim_days,
                   t_max,
                   t_min,
                   humidty,
                   # surface_area,
                   dd_effect,
                   max_survival,
                   loci_cov_matrix,
                   ldt,
                   mu,
                   sigma_dd) {
    # if (sim_days == 25) browser()
    # browser()
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
      # max_temp <- t_max[sim_days, i]
      # min_temp <- t_min[sim_days, i]
      daily_temp <- (t_max+t_min)/2
      
      delay <- sim_delays(n.fem, daily_temp)
      batch_sizes <- sim_batch_sizes(n.fem)
      
      #### Oviposition conditions. This simulates the intervals between when the 
      #### female becomes gravid and oviposition
      cond1 <- as.numeric(fem$next_oviposition >= delay & fem$parity1 == 0 & fem$gravid == 1)
      cond2 <- as.numeric(fem$next_oviposition >= delay & fem$parity1 == 1 & fem$parity2 == 0 & fem$gravid == 1)
      cond3 <- as.numeric(fem$next_oviposition >= delay & fem$parity2 == 1 & fem$parity3 == 0 & fem$gravid == 1)


      
      #Effect of drive on fitness {0 = wild type (w) 1 = drive (d), 2 = non-functional resistance (r2)} 
      
      # wt_homozygous  <- (fem$autosome1 == 0) & (fem$autosome2 == 0)
      # drive_homozygous <- (fem$autosome1 == 1) & (fem$autosome2 == 1)  # sterile
      # res_homozygous <- fem$autosome1 == 2 & fem$autosome2 == 2   # sterile
      # drive_wt <- ((fem$autosome1 == 0) & (fem$autosome2 == 1)) | ((fem$autosome1 == 1) & (fem$autosome2 == 0))
      # drive_res <- ((fem$autosome1 == 1) & (fem$autosome2 == 2)) | ((fem$autosome1 == 2) & (fem$autosome2 == 1))   # sterile
      # res_wt <- ((fem$autosome1 == 0) & (fem$autosome2 == 2)) | ((fem$autosome1 == 2) & (fem$autosome2 == 0))
      # disrupted_loci <- drive_homozygous + res_homozygous + drive_res
      # 

      
      wt_homozygous  <- (fem$autosome1 == 0) & (fem$autosome2 == 0)
      wt2_homozygous  <- (fem$autosome1 == 1) & (fem$autosome2 == 1)
      drive_homozygous <- (fem$autosome1 == 2) & (fem$autosome2 == 2)  # sterile
      res_homozygous <- fem$autosome1 == 3 & fem$autosome2 == 3   # sterile
      drive_wt <- ((fem$autosome1 == 0) & (fem$autosome2 == 2)) | ((fem$autosome1 == 2) & (fem$autosome2 == 0))
      drive_wt2 <- ((fem$autosome1 == 1) & (fem$autosome2 == 2)) | ((fem$autosome1 == 2) & (fem$autosome2 == 1))
      drive_res <- ((fem$autosome1 == 2) & (fem$autosome2 == 3)) | ((fem$autosome1 == 3) & (fem$autosome2 == 2))   # sterile
      res_wt <- ((fem$autosome1 == 0) & (fem$autosome2 == 3)) | ((fem$autosome1 == 3) & (fem$autosome2 == 0))
      res_wt2 <- ((fem$autosome1 == 1) & (fem$autosome2 == 3)) | ((fem$autosome1 == 3) & (fem$autosome2 == 1))
      disrupted_loci <- drive_homozygous + res_homozygous + drive_res
      
      
     # oviposition (with effect of deleterious allele on fitness: sterility)
      
      if (drive_type == "homing" & drive_effect == "sterility") {
        sterile_loci <- (disrupted_loci > 0)
        sterility <- as.numeric(rowSums(sterile_loci) == 0)
        exp_offspring1 <- cond1 * fem$gravid * batch_sizes * sterility
        exp_offspring2 <- cond2 * fem$gravid * batch_sizes * sterility
        exp_offspring3 <- cond3 * fem$gravid * batch_sizes * sterility
        exp_offspring <- exp_offspring1 + exp_offspring2 + exp_offspring3
        
      } else if (drive_type == "homing" & drive_effect == "x_shredding") {
        exp_offspring1 <- cond1 * fem$gravid * batch_sizes
        exp_offspring2 <- cond2 * fem$gravid * batch_sizes
        exp_offspring3 <- cond3 * fem$gravid * batch_sizes
        exp_offspring <- exp_offspring1 + exp_offspring2 + exp_offspring3

      # } else if (drive_type == "toxin_antidote") {
        #   sterile_loci <- (disrupted_loci > 0)
        #   sterility <- as.integer(rowSums(sterile_loci) != ncol(sterile_loci))
        #   exp_offspring1 <- cond1 * fem$gravid * batch_sizes * sterility
        #   exp_offspring2 <- cond2 * fem$gravid * batch_sizes * sterility
        #   exp_offspring3 <- cond3 * fem$gravid * batch_sizes * sterility
        #   exp_offspring <- exp_offspring1 + exp_offspring2 + exp_offspring3
      }
        else {
        exp_offspring1 <- cond1 * fem$gravid * batch_sizes
        exp_offspring2 <- cond2 * fem$gravid * batch_sizes
        exp_offspring3 <- cond3 * fem$gravid * batch_sizes
        exp_offspring <- exp_offspring1 + exp_offspring2 + exp_offspring3
      }
      
      #Update the reproduction memory for all adult females 
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
    
    #Changed NA to "zeros" in case the operation produced NAs from none mated 
    # individuals that are still part of the female population 
    exp_offspring <- replace(exp_offspring, is.na(exp_offspring), 0) 
    
    
    # Offspring generation: 
    
    # n_offspring <- rpois(n.fem, exp_offspring)  
    # No need to redraw from poison since batch sizes has already been calculated using a mean value
    
    total_offspring <- sum(exp_offspring, na.rm = TRUE)
      
    
      if (total_offspring > 0){  
      # Genetic makeup of gametes: replicate these according to the number `exp_offspring` produced per female/male pair
        if (drive_type == "homing" & drive_effect == "sterility"){
          fem_germline <- fem[rep(1:n.fem, exp_offspring), c("autosome1", "autosome2")]
          male_germline <- fem[rep(1:n.fem, exp_offspring), c("male_autosome1", "male_autosome2")]
          fem_germline <- home_drive_conv(fem_germline, chrom1 = "autosome1", chrom2 = "autosome2", prob1, prob2)
          male_germline <- home_drive_conv(male_germline, chrom1 = "male_autosome1", chrom2 = "male_autosome2", prob1, prob2)
          
          fem_allo <- fem[rep(1:n.fem, exp_offspring), c("allo1", "allo2")]
          male_allo <- fem[rep(1:n.fem, exp_offspring), c("male_allo1", "male_allo2")]
        } else if (drive_type == "homing" & drive_effect == "x_shredding"){
          fem_germline <- fem[rep(1:n.fem, exp_offspring), c("autosome1", "autosome2")]
          male_germline <- fem[rep(1:n.fem, exp_offspring), c("male_autosome1", "male_autosome2")]
          fem_germline <- shred_drive_conv(fem_germline, chrom1 = "autosome1", chrom2 = "autosome2", homing_prob)
          male_germline <- shred_drive_conv(male_germline, chrom1 = "male_autosome1", chrom2 = "male_autosome2", homing_prob)

          fem_allo <- fem[rep(1:n.fem, exp_offspring), c("allo1", "allo2")]
          male_allo <- fem[rep(1:n.fem, exp_offspring), c("male_allo1", "male_allo2")]
          
         
          male_allo <- shred_x (male_germline, chrom1 = "male_autosome1", 
                                chrom2 = "male_autosome2", male_allo, 
                                allosome1 = "male_allo1", 
                                allosome2 = "male_allo2", shred_prob)
          
        # } else if (drive_type == "toxin_antidote"){
        #   fem_germline <- fem[rep(1:n.fem, exp_offspring), c("autosome1", "autosome2")]
        #   male_germline <- fem[rep(1:n.fem, exp_offspring), c("male_autosome1", "male_autosome2")]
        #   
        #   fem_allo <- fem[rep(1:n.fem, exp_offspring), c("allo1", "allo2")]
        #   male_allo <- fem[rep(1:n.fem, exp_offspring), c("male_allo1", "male_allo2")]
        #   

        } else {   # this keeps it Mendelian (without recombination)
          fem_germline <- fem[rep(1:n.fem, exp_offspring), c("autosome1", "autosome2")]
          male_germline <- fem[rep(1:n.fem, exp_offspring), c("male_autosome1", "male_autosome2")]
          
          fem_allo <- fem[rep(1:n.fem, exp_offspring), c("allo1", "allo2")]
          male_allo <- fem[rep(1:n.fem, exp_offspring), c("male_allo1", "male_allo2")]
        }
      
        
      
        #check for errors...
        num_loci <- ncol(fem_germline$autosome1)
        stopifnot(num_loci == n_loci)
        
        
        
        # inheritance with or without recombination 
        
        if (recomb) {
        
          cov_matrix <- place_loci_mat(n_loci, genome.size = 1, var = 1, decay)

          # selects which combination of loci based on their proximity are passed to the offspring  
          which_autosome <- which_allele_fn(total_offspring, n_loci, cov_matrix) # female gametes
          which_autosome_mate <- which_allele_fn(total_offspring, n_loci, cov_matrix) # male gametes
        } else {
          
          # this selects which autosome arm is passed to the offspring i.e. no recombination
          which_autosome  <- matrix(rep(rbinom(total_offspring, 1, 0.5), n_loci),
                                    nrow = total_offspring, ncol = n_loci)
          which_autosome_mate <- matrix(rep(rbinom(total_offspring, 1, 0.5), n_loci),
                                        nrow = total_offspring, ncol = n_loci)
        }
 
       # this selects which allosome (sex determining chromosome) arm is passed to the offspring
        which_allosome  <- matrix(rbinom(total_offspring, 1, 0.5), 1,
                                  nrow = total_offspring, ncol = 1)
        which_allosome_mate <- matrix(rbinom(total_offspring, 1, 0.5), 1,
                                      nrow = total_offspring, ncol = 1)
        
        
        
      # Determination of offspring features
        offspring <- tibble(
          stage = "egg",
          allo1 = ifelse(which_allosome,
                         fem_allo$allo1,
                         fem_allo$allo2),           
          allo2 = ifelse(which_allosome_mate,
                         male_allo$male_allo1,
                         male_allo$male_allo2),     
          sex = if_else(allo1[,1] == "X" & allo2[,1] == "X", "female", "male"),
          autosome1 = ifelse(which_autosome,
                             fem_germline$autosome1,
                             fem_germline$autosome2), 
          autosome2 =  ifelse(which_autosome_mate,
                              male_germline$male_autosome1,
                              male_germline$male_autosome2),
          male_autosome1 = matrix(NA, nrow = total_offspring, ncol = n_loci),
          male_autosome2 = matrix(NA, nrow = total_offspring, ncol = n_loci),
          male_allo1 = matrix(NA_character_, nrow = total_offspring, ncol = 1),
          male_allo2 = matrix(NA_character_, nrow = total_offspring, ncol = 1),
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
      
        
        # # effect of drive on individual fitness: lethal effect
        # 
        # if (drive_type == "homing" & lethal_effect){
        #   drive_homozygous <- (offspring$autosome1 == 1) & (offspring$autosome2 == 1)  # sterile
        #   res_homozygous <- offspring$autosome1 == 2 & offspring$autosome2 == 2   # sterile
        #   drive_res <- ((offspring$autosome1 == 1) & (offspring$autosome2 == 2)) | ((offspring$autosome1 == 2) & (offspring$autosome2 == 1))   # sterile
        #   disrupted_loci <- drive_homozygous + res_homozygous + drive_res
        # 
        #   homozygous_lethal <- (disrupted_loci > 0)
        #   any_lethal <- rowSums(homozygous_lethal) > 0
        #   offspring <- filter(offspring, !any_lethal)
        #   #pop <- pop[pop[!any_lethal], ]
        # }


        # Update pop with offspring & fem population
        pop <- pop[!(pop$sex == "female" & pop$stage == "adult"), ]
        pop <- bind_rows(pop, offspring, fem)
      } else {
        # Update pop with females only if no oviposition happened 
        pop <- pop[!(pop$sex == "female" & pop$stage == "adult"), ]
        pop <- bind_rows(pop, fem)
      }
      
    # stage development using growth-degree day accumulation
 
    egg_gdd_accumulated <- cal_dd (t_max, t_min, ldt["egg"])
    larva_gdd_accumulated <- cal_dd (t_max, t_min, ldt["larva"])
    pupa_gdd_accumulated <- cal_dd (t_max, t_min, ldt["pupa"])

    
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
   # population density: aquatic stages (eggs, larvae, and pupae) 
   aq_stage_count <- sum((pop$stage == "egg") + (pop$stage == "larva") + (pop$stage == "pupa"))
   
   # aquatic_stage_density <- aq_stage_count/surface_area
   
   # daily temp
   daily_temp <- (t_max+t_min)/2
   
   # daily humidity  
   # daily_humidity <- humidity

   # #Aquatic stage density-dependent survival using daily mortality hazard and developmental rate (see Golding et al., unpublished data)  
   #   pop <- pop |> mutate(
   #     alive = case_when(
   #       stage == "egg" ~ rbinom(n(), 1, das_temp_dens_As(daily_temp, aquatic_stage_density)),
   #       stage == "larva" ~ rbinom(n(), 1, das_temp_dens_As(daily_temp, aquatic_stage_density)),
   #       stage == "pupa" ~ rbinom(n(), 1, das_temp_dens_As(daily_temp, aquatic_stage_density)),
   #       stage == "adult" ~ rbinom(n(), 1, ds_temp_humid_As(daily_temp, daily_humidity, species = "An. stephensi")),
   #       TRUE ~ NA_integer_
   #     ),
   #     alive = alive == 1
   # )

  
   #Aquatic stage density-dependent survival using Beverton-Holt survival model
   
   pop <- pop |> mutate(
     alive = case_when(
       stage == "egg" ~ rbinom(n(), 1, b_holt_survival(aq_stage_count, max_survival, dd_effect)),
       stage == "larva" ~ rbinom(n(), 1, b_holt_survival(aq_stage_count, max_survival, dd_effect)),
       stage == "pupa" ~ rbinom(n(), 1, b_holt_survival(aq_stage_count, max_survival, dd_effect)),
       stage == "adult" ~ rbinom(n(), 1, max_survival),
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
                       # stages,
                       prob_wildtype1, 
                       prob_wildtype2,
                       per_release,
                       release_day,
                       bloodmeal_prob, 
                       beta,
                       decay,
                       recomb,
                       # lethal_effect,
                       drive_effect,
                       drive_type,
                       prob1,
                       prob2,
                       homing_prob,
                       shred_prob,
                       sim_days,
                       dispersal_type,
                       t_max,
                       t_min,
                       humidty,
                       # surface_area,
                       dd_effect,
                       max_survival,
                       ldt,
                       mu,
                       sigma_dd) {
  #initialise population
  pop <- ini_pop(patches, 
                 n_per_patch, 
                 coords, 
                 n_loci, 
                 # stages,
                 prob_wildtype1, 
                 prob_wildtype2) 

  
  patch_sizes <- list()
  genetic_data <- list()
  
  
  for (day in 1:sim_days) {
   # if (day == 200) browser()
    if (drive_type == "homing" & day == release_day){
      drive_pop <- introduce_drive(initial_pop = nrow(pop[[1]]), per_release, n_loci)
      pop[[1]] <- bind_rows(pop[[1]], drive_pop)
    }
    
    # Growth with reproduction
    pop <- growth(pop_patches = pop,
                  bloodmeal_prob, 
                  n_loci,
                  beta,
                  decay,
                  recomb,
                  # lethal_effect,
                  drive_effect,
                  drive_type,
                  prob1,
                  prob2,
                  homing_prob,
                  shred_prob,
                  sim_days = day,
                  t_max,
                  t_min,
                  humidty,
                  # surface_area,
                  dd_effect,
                  max_survival,
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
      males = sapply(pop, function(p) sum(p$sex == "male", na.rm = TRUE)),
      females = sapply(pop, function(p) sum(p$sex == "female", na.rm = TRUE)),
      pop_size = sapply(pop, nrow),
      patch_occupied = sum(pop_size > 0),
      unoccupied = length(patch) - patch_occupied
    )
    
    
    
    genetic_data[[day]] <- lapply(seq_along(pop), function(patch_id) {
      # browser()
      patch_pop <- pop[[patch_id]]
      autosome1   <- patch_pop$autosome1   
      autosome2   <- patch_pop$autosome2   
      n_ind   <- nrow(autosome1)
      n_loci  <- ncol(autosome1)
      
      total_alleles <- 2 * n_ind
      wildtype1 <- colSums(autosome1 == 0) + colSums(autosome2 == 0)
      wildtype2 <- colSums(autosome1 == 1) + colSums(autosome2 == 1)
      
      drive <- colSums(autosome1 == 2) + colSums(autosome2 == 2)
      resistance <- colSums(autosome1 == 3) + colSums(autosome2 == 3)
      
      freq_w1 <- ifelse(wildtype1 > 0, wildtype1 / total_alleles, 0)
      freq_w2 <- ifelse(wildtype2 > 0, wildtype2 / total_alleles, 0)
      freq_d <- ifelse(drive > 0, drive / total_alleles, 0)
      freq_r <- ifelse(resistance > 0, resistance / total_alleles, 0)
      
      tibble(
        patch = patch_id,
        time_step  = day,
        locus = 1:n_loci,
        w1 = freq_w1,
        w2 = freq_w2,
        d = freq_d,
        r = freq_r
      )
    })
    
    patch_sizes_df <- bind_rows(patch_sizes)
    genetic_data_df <- bind_rows(genetic_data)
    cat("Day", day, "Completed \n")
    
  }
  
  # Return the collected data
  list(
    pop_sizes = patch_sizes_df,
    genetic_df = genetic_data_df,
    final_pop = pop
  )
}


