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
      gdd_accumulated = 0,
      parity1 = 0,
      parity2 = 0,
      alive = TRUE
    )
    if (length(n_per_patch) != patches) warning("Initial patch population does not equal specified number of patches")
  }
  
  return(patches_pop)
}



# Loci selection matrix: function to place loci at random on the genome (of size = 1)
# also takes exponential decay and variance to produce variance-covariance matrix

place_loci_mat <- function(loci, genome.size = 1, var = 1, decay){
  loci_positions <- sort(runif(loci, max = genome.size))
  loci_dist_matrix <- as.matrix(dist(loci_positions))^2 
  loci_cov_matrix <- var*exp(-decay*loci_dist_matrix)
  return(loci_cov_matrix)
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
                   sim_days,
                   t_max,
                   t_min,
                   sigma,
                   loci_cov_matrix,
                   gdd_required,
                   ldt) {
  
#if (sim_days == 50) browser()
  updated_pop_patches <- list()
  
  for (i in seq_along(pop_patches)) {
    pop <- pop_patches[[i]]  
    
    male <- pop |> filter(sex == 0, stage == "adult")    # All males
    fem <- pop |> filter(sex == 1, stage == "adult")     # All females
    n.fem <- nrow(fem)
    n.male <- nrow(male)

    if (n.fem > 0 && n.male > 0){
      # select a mate for each female
      selected_male_index <- sample(x = n.male, size = n.fem, replace = TRUE)
      selected_male <- male[selected_male_index, ]
      
      
      # Bernoulli trial for mating and feeding (1 if mated/bloodfed, 0 if not)
      realised_mated <- rbinom(n = n.fem, 1, prob = (n.male/(beta + n.male)))#  (nrow(male)/(beta + nrow(male)))) is the mating probability which increases as male population increases (North and Godfray; Malar J (2018) 17:140) 
      realised_bloodmeal <- rbinom(n = n.fem, 1, bloodmeal_prob)   
      gravid <- as.integer(realised_mated*realised_bloodmeal == 1)
      
      # effect of load on fecundity: to turn this effect off, set fecundity_effect = 0 in function call
      # If bloodfed, calculate expected offspring for this female
      # Additive effect as the number of homozygous lethal loci increases fecudity 
      # decreases exponentially as given below:
      
      homo_loci <- apply(fem$allele1 + fem$allele2 == 2, 1, sum) # number of homozygous loci for each female
      exp_offspring <- gravid * fecundity * exp(-fecundity_effect * homo_loci)
      
      # # complete sterility effect
      # homo_loci <- 1 - (as.integer(apply(fem$allele1 + fem$allele2 == 2, 1, any))) # if any loci is homozygous
      # exp_offspring <- gravid * fecundity * homo_loci

      
    # oviposition based on degree-day accumulation
 
        ##### underway###
      
    
      # Draw the actual number of offspring from a Poisson distribution
      n_offspring <- rpois(n = n.fem, exp_offspring) 
      
      
      }  else {
        # If not, set offspring count to 0
        n_offspring <- rep(0, n.fem)
      }
  
      
      total_offspring <- sum(n_offspring)
      
      if (total_offspring > 0){
      # Replicate the parents features `n_offspring` times for each offspring, collect only genetic information
        fem_germline <- fem[rep(1:n.fem, n_offspring), ] |> select(contains("allele"))
        male_germline <- selected_male[rep(1:n.fem, n_offspring), ] |> select(contains("allele"))
  
        # Genetic inheritance
        num_loci <- ncol(fem_germline$allele1)
        stopifnot(num_loci == n_loci)
      
        # random selection of allele, with linkage 
        which_allele_fn <- function(n_offspring, num_loci, loci_cov_matrix){
          epsilon <- MASS::mvrnorm(n_offspring, rep(0, num_loci), Sigma = loci_cov_matrix)
          selection_prob <- plogis(epsilon)
          matrix(rbinom(n_offspring * num_loci, 1, selection_prob) == 1,
               nrow = n_offspring,
               ncol = num_loci)
        }
      
        which_allele_female <- which_allele_fn(total_offspring, num_loci, loci_cov_matrix) #female gametes
        which_allele_male <- which_allele_fn(total_offspring, num_loci, loci_cov_matrix) # male gametes
      
        #  Determination of offspring features
        offspring <- tibble(
          sex = rbinom(total_offspring, 1, 0.5),
          stage = "egg",
          allele1 = ifelse(which_allele_female,
                         fem_germline$allele1,
                         fem_germline$allele2),
          allele2 = ifelse(which_allele_male,
                         male_germline$allele1,
                         male_germline$allele2),
          gdd_accumulated = 0,
          parity1 = 0,
          parity2 = 0,
          alive = TRUE
        )
      
        # Add offspring to main population
        pop <- bind_rows(pop, offspring)
      }
      
      # Genetic load: lethal effect
      
      if (lethal_effect){
        homozygous_lethal <- (pop$allele1 == 1) & (pop$allele2 == 1)
        any_homozygous <- rowSums(homozygous_lethal) > 0
        pop <- filter(pop, !any_homozygous)
      }
      
      
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
  
      
    
    # Temperature variable  
    
      max_temp <- t_max[i]
      min_temp <- t_min[i]
      #daily_temp <- (max_temp+min_temp)/2
      

    # Temperature-adjusted survival for larval and adult population: NOTE: set 
    #  SD in "temp" to "0" to turn off temperature variation onn survival
    
    # temp_effect <- exp(-((daily_temp - 25)^2) / (2 * sigma^2)) # Gaussian process (modified from Beck-Johnson et al., 2013). 
    # temp_adjusted_survival <- daily_survival[c("adult","larva")] * temp_effect
    
    
    # Density-dependent survival for larval stage
    larva_count <- sum(pop$stage == "larva")
    density_dependence <- 1/(1 + (alpha*larva_count))
    #density_dependent_survival <- temp_adjusted_survival["larva"] * density_dependence
    density_dependent_survival <- daily_survival["larva"] * density_dependence
    
    pop <- pop |> mutate(
        alive = case_when(
          stage == "egg" ~ rbinom(n(), 1, daily_survival["egg"]),
          stage == "larva" ~ rbinom(n(), 1, density_dependent_survival),
          stage == "pupa" ~ rbinom(n(), 1, daily_survival["pupa"]),
          stage == "adult" ~ rbinom(n(), 1, daily_survival["adult"]),
        ),
        alive = alive == 1
     )
    
    # function to calculate growth degree-days
    gdd_cal <- function(daily_max_temp, daily_min_temp, T_base) {
      pmax(0, (daily_max_temp+daily_min_temp)/2 - T_base)
    }
  
    daily_gdd_accumulated <- gdd_cal (max_temp, min_temp, ldt)
    pop <- pop |> mutate(gdd_accumulated = gdd_accumulated + daily_gdd_accumulated)    
    
    
# stage transition probability per time step based on degree day accumlated

 pop <- pop |>
      mutate(
        transition_egg_larva  = stage == "egg"  & gdd_accumulated >= gdd_required["egg_larva"]  & rbinom(n(), 1, daily_transition["egg_larva"]) == 1,
        transition_larva_pupa = stage == "larva" & gdd_accumulated >= gdd_required["larva_pupa"] & rbinom(n(), 1, daily_transition["larva_pupa"]) == 1,
        transition_pupa_adult = stage == "pupa"  & gdd_accumulated >= gdd_required["pupa_adult"] & rbinom(n(), 1, daily_transition["pupa_adult"]) == 1
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
                       sim_days,
                       dispersal_matrix,
                       t_max,
                       t_min,
                       sigma,
                       gdd_required,
                       ldt) {
  pop <- ini_pop(patches, n_per_patch, coords, n_loci, init_frequency)
  l.cov.mat <- place_loci_mat(n_loci, genome.size = 1, var = 1, decay)
  
  patch_sizes <- list()
  stage_distributions <- list()
  for (day in 1:sim_days) {
    #if (day == 9) browser()
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
                  sim_days = day,
                  t_max,
                  t_min,
                  sigma,
                  loci_cov_matrix = l.cov.mat,
                  gdd_required,
                  ldt)
    
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


