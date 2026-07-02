

# function to load required packages 

load_libraries <- function(pack) {
  for (p in pack) {
    if (!requireNamespace(p, quietly = TRUE)) {
      install.packages(p)
    }
    library(p, character.only = TRUE)
  }
}


#initialise individuals in patches 
create_n_per_patch <- function(patches, n_individual) {
  
  if (patches < 1) {
    stop("Number of patches must be at least 1.")
  }
  if (n_individual < 0) {
    stop("Carrying capacity must not be < 0")
  }
  n_per_patch <- rep(0, patches)
  n_per_patch[1] <- n_individual
  return(n_per_patch)
}





# # function to make autosome
# 
# make_autosome <- function(individual, prefix, n_loci) {
#  
#    matrix(
#   c(paste0(rep(prefix, individual * (n_loci)), rep(1:(n_loci), each = individual))),
#   nrow = individual, ncol = n_loci
# )
# 
# }
# 
# #make_autosome(5, "A", 5)



# function to make allosome (sex chromosome)

# make_allosome <- function(sex_alleles, individual) {
#   
#   cbind(
#     sample(sex_alleles, size = individual, replace = TRUE),
#     rep(0L, individual)
#   )
# }
#   
# 
# make_allosome <- function(sex_alleles, individual) {
#   sample(sex_alleles, size = individual, replace = TRUE)
#   }



make_allosome <- function(sex_alleles, individual) {
  matrix(sample(sex_alleles, size = individual, replace = TRUE), 
         nrow =individual, ncol = 1)
}


# correlated allele selection (recombination) 

# Loci selection matrix: function to randomly assign position to all loci on the genome (of size = 1)
# also takes exponential decay and variance to produce variance-covariance matrix

place_loci_mat <- function(loci, genome.size = 1, var = 1, decay){
  loci_positions <- sort((runif(loci, max = genome.size)))
  loci_dist_matrix <- as.matrix(dist(loci_positions))^2 
  loci_cov_matrix <- var*exp(-decay*loci_dist_matrix)
  return(loci_cov_matrix)
}


# function to generate random multivariate normal effect values and transform them to probabilities
which_allele_fn <- function(exp_offspring, num_loci, loci_cov_matrix){
  epsilon <- MASS::mvrnorm(exp_offspring, rep(0, num_loci), Sigma = loci_cov_matrix)
  selection_prob <- plogis(epsilon)
  matrix(rbinom(exp_offspring * num_loci, 1, selection_prob) == 1,
         nrow = exp_offspring,
         ncol = num_loci)
}


# alternative  function maybe faster (computational speed)
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



# growth degree day estimation (Abbasi et al., Environmental Entomology, 2023, Vol. 52, No. 6)

# function to calculate growth degree-days accumulation 

cal_dd <- function(daily_max_temp, daily_min_temp, T_base) {
  max(0, (daily_max_temp+daily_min_temp)/2 - T_base)
}


# function to estimate the probability of transition based on mean growth degree 
# days at different percentile using inverse cumulative distribution function* 

erfinv <- function (x) qnorm((1 + x)/2)/sqrt(2)

sigma_etimate <- function(x, mu, p) {
  (x - mu)/(sqrt(2) * erfinv(2*p-1))
}

prob_trans <- function(dd, mu, sigma) {
  pnorm(dd, mu, sigma)
}


# # estimated survival based on temperature and population density (aquatic stage),  
# # and temperature and humidity (adult stage) using daily mortality hazard, and 
# # developmental rate (see Golding et al., unpublished data)
# 
# ensure_positive <- function(x) {
#   x * as.numeric(x > 0)
# }
# 
# # reload lifehistory functions from saved objects (RDS file) to used for survival. 
# # Adapted from Golding et al., unpublished)
# 
# rehydrate_lifehistory_function <- function(path_to_object) {
#   object <- readRDS(path_to_object)
#   do.call(`function`,
#           list(object$arguments,
#                body(object$dummy_function)))
# }
# 
# 
# aquatic_stage <- "C:/Users/JOhiolei/OneDrive - The Kids Research Institute Australia/Documents/An-stephensi-IBM_Proj/R/das_temp_dens_As.RDS"
# adult_stage <- "C:/Users/JOhiolei/OneDrive - The Kids Research Institute Australia/Documents/An-stephensi-IBM_Proj/R/ds_temp_humid.RDS"
# 
# das_temp_dens_As <- rehydrate_lifehistory_function(aquatic_stage)
# ds_temp_humid_As <- rehydrate_lifehistory_function(adult_stage)





# Density-dependent survival for the aquatic stage 
b_holt_survival <- function(density, 
                              max_survival = 0.85,
                              dd_effect = 0.01) {
  # Or Zimmermann et al., 2021; ICES Journal of Marine Science (2021), 78(6) 2193 2203. doi:10.1093/icesjms/fsaa246
   max_survival / (1 + dd_effect * density) # density = number of aquatic stage individuals
}





# function to simulate oviposition frequency and batch sizes. mean eggs per female
# per day (EFD) and mean temperature were estimated from Villena et al., https://doi.org/10.1002/ecy.3685

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


# We simulate the expected batch sizes based on Suleman, 1990 https://doi.org/10.1093/jmedent/27.5.819
# to match the mean and SD but modelled as negative binomial
sim_batch_sizes <- function(n) {
  rnbinom(n, mu = 96.8, size = 1 / 0.16)
}

# calculate the frequency of oviposition or expected delay between batches,
# given the expected batch size
# 
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

expected_egg_laying_delay <- function(temp, expected_batch_size = 96.8) {
  expected_batch_size / egg_laying_rate(temp)
}

sim_delays <- function(n, temp) {
  expected_delay <- expected_egg_laying_delay(temp)
  delays_continuous <- rlnorm_mean_var(n,
                                       expected_delay,
                                       sd = 0.5)
  delays <- pmax(1, round(delays_continuous))
  delays
}



#### negative exponential dispersal kernel  

metapop <- function(coords, lambda, disp_prob) {
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
  dispersal_matrix <- disp_prob * rel_dispersal_matrix +
    (1 - disp_prob) * diag(nrow(dispersal_kernel))
  
  return(dispersal_matrix)
}


# adjacency matrix 

step_stone <- function(n_patches, disp_prob) {
  
  matrix_landscape <- matrix(0, n_patches, n_patches)
  adjacency <- abs(row(matrix_landscape) - col(matrix_landscape)) == 1
  adjacency[] <- as.numeric(adjacency)

  # make these rows sum to 1 to get probability of moving to other patch
  # *if* they left. This dispersal matrix gives the probability of the vector
  # vector moving between patches
  rel_dispersal_matrix <- sweep(adjacency, 1,
                                rowSums(adjacency), FUN = "/")

  # normalise these to have the overall probability of dispersing to that patch,
  # and add back the probability of remaining
  dispersal_matrix <- disp_prob * rel_dispersal_matrix +
    (1 - disp_prob) * diag(nrow(adjacency))

  return(dispersal_matrix)
}






#Homing gene drive function (conversion mechanism)

home_drive_conv <- function(parent, chrom1, chrom2, prob1, prob2) {
  # browser()
  
  loci1 <- parent[[chrom1]] 
  loci2 <- parent[[chrom2]] 
  
  # if (any(is.na(loci1)) | any(loci2)) {
  #     warning("NA detected in allele input!")
  # }
  
  drive_wt <- ((loci1 == 0) & (loci2 == 2)) | ((loci1 == 2) & (loci2 == 0))
  
  drive_wt <- 1 * drive_wt  # convert logical to numeric
  
  
  #cleavage
  cleavage  <- matrix(rbinom(nrow(loci1), 1, prob1), ncol(loci1), # drive cleavage at each locus
                      nrow = nrow(loci1), ncol = ncol(loci1))
  # homing
  homing  <- matrix(rbinom(nrow(loci1), 1, prob2), ncol(loci1), # drive conversion at each locus
                    nrow = nrow(loci1), ncol = ncol(loci1)) 
  conv_event <- homing*cleavage # conversion event?
  
  loci1[loci1 == 0 & conv_event == 1 & drive_wt == 1] <- 2 # successful conversions
  loci2[loci2 == 0 & conv_event == 1 & drive_wt == 1] <- 2
  
  # # Resistance development if homing fails (0 to 2)
  # resist_dev  <- matrix(rbinom(nrow(loci1), 1, 1-prob2), ncol(loci1), # drive conversion at each locus
  #                       nrow = nrow(loci1), ncol = ncol(loci1))
  # res_event <- resist_dev *cleavage # resistance development through NHEJ
  # loci1[loci1 == 0 &  res_event == 1 & drive_wt == 1] <- 3  
  # loci2[loci2 == 0 & res_event == 1 & drive_wt == 1] <- 3
  
  loci1[loci1 == 0 &  conv_event == 0 & drive_wt == 1] <- 3  
  loci2[loci2 == 0 & conv_event == 0 & drive_wt == 1] <- 3
  
  parent[[chrom1]]   <- loci1
  parent[[chrom2]]   <- loci2
  
  return(parent)
}


#sex-shredder homing function
shred_drive_conv <- function(parent, chrom1, chrom2, homing_prob) {
  # browser()

  loci1 <- parent[[chrom1]]
  loci2 <- parent[[chrom2]]

  # if (any(is.na(loci1)) | any(loci2)) {
  #     warning("NA detected in allele input!")
  # }

  drive_wt <- ((loci1 == 0) & (loci2 == 2)) | ((loci1 == 2) & (loci2 == 0))

  drive_wt <- 1 * drive_wt  # convert logical to numeric

  # homing
  homing  <- matrix(rbinom(nrow(loci1), 1, homing_prob), ncol(loci1), # drive conversion at each locus
                    nrow = nrow(loci1), ncol = ncol(loci1))

  loci1[loci1 == 0 & homing == 1 & drive_wt == 1] <- 2 # successful conversions
  loci2[loci2 == 0 & homing == 1 & drive_wt == 1] <- 2


  parent[[chrom1]]   <- loci1
  parent[[chrom2]]   <- loci2

  return(parent)
}


# this function follows the #sex-shredder homing function to complete the sex bias in the germline
shred_x <- function(parent_auto, chrom1, chrom2, parent_allo, allosome1, allosome2, shred_prob) {

  auto1 <- parent_auto[[chrom1]]
  auto2 <- parent_auto[[chrom2]]

  allo1 <- parent_allo[[allosome1]]
  allo2 <- parent_allo[[allosome2]]

  drive <- (auto1 == 2) & (auto2 == 2) 
  
  drive <- 1 * drive  # convert logical to numeric

  # shredding
  shredding  <- matrix(rbinom(nrow(allo1), 1, shred_prob), ncol(allo1), # shredding
                       nrow = nrow(allo1), ncol = ncol(allo1))

  allo1[allo1 == "X" & allo2 == "Y" & drive == 1 & shredding == 1] <- "Y" #shredding of X

  parent_allo[[allosome1]] = allo1

  return(parent_allo)
}



# #sex-shredder homing function
# shred_drive_conv <- function(parent_auto, auto1, auto2, homing_prob, parent_allo, allo1, allo2, shred_prob) {
#   # browser()
#   
#   loci_1 <- parent_auto[[auto1]] 
#   loci_2 <- parent_auto[[auto2]] 
#   
#   loci_a <- parent_allo[[allo1]]
#   loci_b <- parent_allo[[allo2]]
#   
#   
#   # if (any(is.na(loci1)) | any(loci2)) {
#   #     warning("NA detected in allele input!")
#   # }
#   
#   drive_wt <- ((loci_1 == 0) & (loci_2 == 2)) | ((loci_1 == 2) & (loci_2 == 0))
#   
#   drive_wt <- 1 * drive_wt  # convert logical to numeric
#   
#   # homing
#   homing  <- matrix(rbinom(nrow(loci_1), 1, homing_prob), ncol(loci_1), # drive conversion at each locus
#                     nrow = nrow(loci_1), ncol = ncol(loci_1)) 
#   
#   loci_1[loci_1 == 0 & homing == 1 & drive_wt == 1] <- 2 # successful conversions
#   loci_2[loci_2 == 0 & homing == 1 & drive_wt == 1] <- 2
#   
#   
#   # shredding
#   shredding  <- matrix(rbinom(nrow(loci_a), 1, shred_prob), ncol(loci_a), # shredding
#                        nrow = nrow(loci_a), ncol = ncol(loci_a)) 
#   
#   loci_a[loci_a == "X" & loci_b == "Y" & drive == 1 & homing == 1 & shredding == 1] <- "Y" #shredding of X
#   
#   
#   
#   parent_auto[[auto1]] <- loci_1
#   parent_auto[[auto2]] <- loci_2
#   parent_allo[[allo1]] <- loci_a
#   parent_allo[[allo2]] <- loci_b
#   
#   return(list(parent_auto, parent_allo))
# }



