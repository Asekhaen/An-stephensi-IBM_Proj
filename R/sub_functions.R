# These functions estimates/calculates parameters that were used in the core functions



# Loci selection matrix: function to place loci at random on the genome (of size = 1)
# also takes exponential decay and variance to produce variance-covariance matrix

place_loci_mat <- function(loci, genome.size = 1, var = 1, decay){
  loci_positions <- sort((runif(loci, max = genome.size)))
  loci_dist_matrix <- as.matrix(dist(loci_positions))^2 
  loci_cov_matrix <- var*exp(-decay*loci_dist_matrix)
  return(loci_cov_matrix)
}


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


# estimated survival based on temperature and population density (aquatic stage),  
# and temperature and humidity (adult stage) using daily mortality hazard, and 
# developmental rate (see Golding et al., unpublished data)

ensure_positive <- function(x) {
  x * as.numeric(x > 0)
}

# reload lifehistory functions from saved objects (RDS file) to used for survival. 
# Adapted from Golding et al., unpublished)

rehydrate_lifehistory_function <- function(path_to_object) {
  object <- readRDS(path_to_object)
  do.call(`function`,
          list(object$arguments,
               body(object$dummy_function)))
}


aquatic_stage <- "C:/Users/22181916/Documents/Curtin-PhD/R_and_IBM/An-stephensi-IBM_Proj/R/das_temp_dens_As.RDS"
adult_stage <- "C:/Users/22181916/Documents/Curtin-PhD/R_and_IBM/An-stephensi-IBM_Proj/R/ds_temp_humid.RDS"

das_temp_dens_As <- rehydrate_lifehistory_function(aquatic_stage)
ds_temp_humid_As <- rehydrate_lifehistory_function(adult_stage)


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



#### create dispersal matrix, called in the core dispersal function (metapopulation)  

make_dispersal_matrix <- function(coords, lambda, dispersal_prop) {
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
  dispersal_matrix <- dispersal_prop * rel_dispersal_matrix +
    (1 - dispersal_prop) * diag(nrow(dispersal_kernel))
  
  return(dispersal_matrix)
}

