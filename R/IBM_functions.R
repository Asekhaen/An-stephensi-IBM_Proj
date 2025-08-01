# Functions for running the stephensi IBM


# create coordinates for the patches/locations 
coords <- as.data.frame(100 * matrix(runif(patches * 2), ncol = 2))
colnames(coords) <- c("x","y")


#### Initialise population ####

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


# Loci selection matrix: function to place loci at random on the genome (of size = 1)
# also takes exponential decay and variance to produce variance-covariance matrix

place_loci_mat <- function(loci, genome.size = 1, var = 1, decay){
  loci_positions <- (runif(loci, max = genome.size))
  loci_dist_matrix <- as.matrix(dist(loci_positions))^2 
  loci_cov_matrix <- var*exp(-decay*loci_dist_matrix)
  return(loci_cov_matrix)
}


# growth degree day estimation (Abbasi et al., Environmental Entomology, 2023, Vol. 52, No. 6)

erfinv <- function (x) qnorm((1 + x)/2)/sqrt(2)
sigma_etimate <- function(x, mu, p) {
  (x - mu)/(sqrt(2) * erfinv(2*p-1))
}

# stage development using degree-days (dd): degree-days, percentage/probability of 
# transitioning (data from Abbasi et al., 2023)   

# Eggs transition
one_per_egg <- sigma_etimate(30.5, 44.5, 0.01)
ten_per_egg <- sigma_etimate(33.7, 44.5, 0.1)
eighty_per_egg <- sigma_etimate(59.6, 44.5, 0.8)
mean_sigma_egg <- mean(c(one_per_egg, ten_per_egg, eighty_per_egg))

#Larva transition
one_per_larva <- sigma_etimate(93.9, 145.5, 0.01)
ten_per_larva <- sigma_etimate(104.1, 145.5, 0.1)
eighty_per_larva <- sigma_etimate(182.3, 145.5, 0.8)
mean_sigma_larva <- mean(c(one_per_larva, ten_per_larva, eighty_per_larva))

#Pupa transition
one_per_pupa <- sigma_etimate(17.7, 29.7, 0.01)
ten_per_pupa <- sigma_etimate(22.3, 29.7, 0.1)
eighty_per_pupa <- sigma_etimate(40.4, 29.7, 0.8)
mean_sigma_pupa <- mean(c(one_per_pupa, ten_per_pupa, eighty_per_pupa))


mu <- c(egg = 44.5, larva = 145.5, pupa = 29.7)
sigma_dd <- c(egg = mean_sigma_egg, 
              larva = mean_sigma_larva, 
              pupa = mean_sigma_pupa)

prob_trans <- function(dd, mu, sigma) {
 pnorm(dd, mu, sigma)
}

# function to calculate growth degree-days 

cal_dd <- function(daily_max_temp, daily_min_temp, T_base) {
  max(0, (daily_max_temp+daily_min_temp)/2 - T_base)
}



# estimated survival based on temperature and population density (aquatic stage),  
# and temperature and humidity (adult stage) using daily mortality hazard, and 
# developmental rate (see Golding et al., unpublished data)

ensure_positive <- function(x) {
  x * as.numeric(x > 0)
}

# reload lifehistory functions from saved objects
rehydrate_lifehistory_function <- function(path_to_object) {
  object <- readRDS(path_to_object)
  do.call(`function`,
          list(object$arguments,
               body(object$dummy_function)))
}

path_aquatic <- "C:/Users/22181916/Documents/Curtin-PhD/R_and_IBM/An-stephensi-IBM_Proj/R/das_temp_dens_As.RDS"
das_temp_dens_As <- rehydrate_lifehistory_function(path_aquatic)


path_adult <- "C:/Users/22181916/Documents/Curtin-PhD/R_and_IBM/An-stephensi-IBM_Proj/R/ds_temp_humid.RDS"
ds_temp_humid_As <- rehydrate_lifehistory_function(path_adult)




# function to simulate oviposition frequency and batch sizes. mean eggs per female
# per day (EFD) and mean temperature were estimated from Villena et al., https://doi.org/10.1002/ecy.3685

egg_laying_rate <- function(temp) {
  peak_temp <- 28
  peak_val <- 26.2
  temp_sd <- 6
  
  unscaled_value <- dnorm(temp,
                          mean = peak_temp,
                          sd = temp_sd)
  normalisation <- dnorm(peak_temp,
                         mean = peak_temp,
                         sd = temp_sd)
  peak_val * unscaled_value / normalisation
  
}

# return the parameters of lognormal with specified mean and variance
lognormal_params <- function(mean, sd) {
  var <- sd ^ 2
  list(
    meanlog = log((mean ^ 2) / sqrt(var + mean ^ 2)),
    sdlog = sqrt(log(1 + var / (mean ^ 2)))
  )
}

# simulate from a lognormal, given the mean and sd of the distribution (not the
# meanlog and sdlog parameters)
rlnorm_mean_var <- function(n, mean, sd) {
  params <- lognormal_params(mean, sd)
  rlnorm(n, params$meanlog, params$sdlog)
}

# simulate delays between egg batches, in days
sim_delays <- function(n, temp) {
  expected_delay <- expected_egg_laying_delay(temp)
  delays_continuous <- rlnorm_mean_var(n,
                                       expected_delay,
                                       sd = 0.5)
  delays <- pmax(1, round(delays_continuous))
  delays
}

# We simulate the expected batch sizes based on Suleman, 1990 https://doi.org/10.1093/jmedent/27.5.819
# to match the mean and SD but modelled as negative binomial
sim_batch_sizes <- function(n) {
  rnbinom(n, mu = 96.8, size = 1 / 0.16)
}

# calculate the frequency of oviposition or expected delay between batches,
# given the expected batch size
expected_egg_laying_delay <- function(temp, expected_batch_size = 96.8) {
  expected_batch_size / egg_laying_rate(temp)
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
                   fecundity_effect,
                   lethal_effect,
                   complete_sterile,
                   sim_days,
                   t_max,
                   t_min,
                   humidty,
                   surface_area,
                   loci_cov_matrix,
                   gdd_required,
                   ldt,
                   mu,
                   sigma_dd) {
 #if (sim_days == 15) browser()
  updated_pop_patches <- list()
  
  for (i in seq_along(pop_patches)) {
    pop <- pop_patches[[i]]  
    
    male <- pop |> filter(sex == 0, stage == "adult")    # All males
    fem <- pop |> filter(sex == 1, stage == "adult")     # All females
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
      
 
      # average daily temperature  
      max_temp <- t_max[i]
      min_temp <- t_min[i]
      daily_temp <- (max_temp+min_temp)/2
      
      delay <- sim_delays(n.fem, daily_temp)
      batch_sizes <- sim_batch_sizes(n.fem)
      
      
      #### Oviposition conditions 

      cond1 <- as.numeric(fem$next_oviposition >= delay & fem$parity1 == 0 & fem$gravid == 1)
      cond2 <- as.numeric(fem$next_oviposition >= delay & fem$parity1 == 1 & fem$parity2 == 0 & fem$gravid == 1)
      cond3 <- as.numeric(fem$next_oviposition >= delay & fem$parity2 == 1 & fem$parity3 == 0 & fem$gravid == 1)

      homo_loci <- rowSums((fem$allele1 + fem$allele2) == 2)        # homozygous loci for each female
      
      
      
      if (complete_sterile) {
        homozygous <- (homo_loci > 0)
        sterile <- as.numeric(!homozygous)
        exp_offspring1 <- cond1 * fem$gravid * batch_sizes * sterile
        exp_offspring2 <- cond2 * fem$gravid * batch_sizes * sterile
        exp_offspring3 <- cond3 * fem$gravid * batch_sizes * sterile
        exp_offspring <- exp_offspring1 + exp_offspring2 + exp_offspring3
      } else {
        exp_offspring1 <- cond1 * fem$gravid * batch_sizes * exp(-fecundity_effect * homo_loci)
        exp_offspring2 <- cond2 * fem$gravid * batch_sizes * exp(-fecundity_effect * homo_loci)
        exp_offspring3 <- cond3 * fem$gravid * batch_sizes * exp(-fecundity_effect * homo_loci)
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
      # If not, set offspring count to 0
      exp_offspring <- rep(0, n.fem)
    }
    
    
    # Offspring generation: Draw the actual number of offspring from a Poisson distribution
    n_offspring <- rpois(n.fem, exp_offspring) # (Sounds logical to use Poisson after estimation based on temperature?)
    total_offspring <- sum(n_offspring)
      
    
      if (total_offspring > 0){  
      # Replicate the parents features `n_offspring` times for each offspring, collect only genetic information

        fem_germline <- fem[rep(1:n.fem, n_offspring), ] |> select(contains("allele"))
        male_germline <- fem[rep(1:n.fem, n_offspring), ] |> select(contains("male"))
         
        
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
        pop <- pop |> filter(!(sex == 1 & stage == "adult"))
        pop <- bind_rows(pop, offspring, fem)
      } else {
        # Update pop with females only
        pop <- pop |> filter(!(sex == 1 & stage == "adult"))
        pop <- bind_rows(pop, fem)
      }
      
    
      # Genetic load: lethal effect
      
      if (lethal_effect){
        homozygous_lethal <- (pop$allele1 == 1) & (pop$allele2 == 1)
        any_homozygous <- rowSums(homozygous_lethal) > 0
        pop <- filter(pop, !any_homozygous)
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
   
   
   #Density dependent survival adjusted to temperature (for aquatic stage) and humidty (for adult stage)
   
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
   

    pop <- filter(pop, alive)
    
    updated_pop_patches[[i]] <- pop
  }
  return(updated_pop_patches)
 }


# Dispersal: the default is the metapopulation network. Otherwise this can be switched to a stepping stone model
# 
# 
# 

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



#### Stepping stone dispersal (one-dimensional discrete space) ####
  

    ss_dispersal <- function(pop_patches, dispersal_frac) {
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
        # can_disperse <- rbinom(nrow(current_patch[current_patch$stage == "adult", ]), 1, dispersal_frac)
        non_adults <- current_patch[current_patch$stage != "adult", ]
        
        # Determine which adults disperse
        ready_to_disperse <- rbinom(nrow(adults), 1, dispersal_frac)
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

  

#### Simulation: Bring them all together ####
simulation <- function(patches,
                       n_per_patch, 
                       coords,
                       n_loci,
                       init_frequency,
                       bloodmeal_prob, 
                       fecundity, 
                       conversion_prob,
                       resistance_prob,
                       daily_survival, 
                       daily_transition,
                       alpha,
                       beta,
                       decay,
                       fecundity_effect,
                       lethal_effect,
                       complete_sterile,
                       sim_days,
                       stepping_stone_model,
                       dispersal_matrix,
                       t_max,
                       t_min,
                       humidty,
                       surface_area,
                       gdd_required,
                       ldt,
                       mu,
                       sigma_dd) {
  pop <- ini_pop(patches, n_per_patch, coords, n_loci, init_frequency)
  l.cov.mat <- place_loci_mat(n_loci, genome.size = 1, var = 1, decay)
  
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
                  fecundity,
                  conversion_prob,
                  resistance_prob,
                  n_loci,
                  daily_survival,
                  daily_transition,
                  alpha,
                  beta,
                  decay,
                  fecundity_effect,
                  lethal_effect,
                  complete_sterile,
                  sim_days = day,
                  t_max,
                  t_min,
                  humidty,
                  surface_area,
                  loci_cov_matrix = l.cov.mat,
                  gdd_required,
                  ldt,
                  mu,
                  sigma_dd)
    
    # Dispersal
    if(stepping_stone_model) {
      pop <- ss_dispersal(pop, dispersal_frac)
    } else {
      pop <- meta_dispersal(pop, dispersal_matrix, check = FALSE)
    }
    
    
    # # Track the average generation time: from egg to first oviposition
    # 
    # generation_time[[day]] <- do.call(rbind, lapply(pop, function(df) {
    #   filter(df, sex == 1, !is.na(first_ovip_day), !is.na(birth_day))
    # }))
    # 
    # generation_time$time_to_first_ovip <- generation_time$first_ovip_day - generation_time$birth_day
    # mean(generation_time$time_to_first_ovip)
    
    
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


