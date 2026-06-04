epsilon <- rnorm(1e4, 0, 1)
p <- plogis(epsilon)
hist(p, xlim = c(0, 1))
y <- rbinom(length(p), size = 1, prob = p)
hist(y)



which_allele_fn <- function(n_offspring, num_loci, loci_cov_matrix){
  epsilon <- MASS::mvrnorm(n_offspring, rep(0, num_loci), Sigma = loci_cov_matrix)
  selection_prob <- plogis(epsilon)
  matrix(rbinom(n_offspring * num_loci, 1, selection_prob) == 1,
         nrow = n_offspring,
         ncol = num_loci)
}

n_loci <- 2
cov <- place_loci_mat(n_loci, decay = 10)
epsilon <- MASS::mvrnorm(1000, rep(0, n_loci), Sigma = cov)
y <- epsilon > 0
# y <- which_allele_fn(1000, n_loci, cov)

cov2cor(cov)[1, 2]

cor(y)[1, 2]


#modified 
which_allele_fn <- function(n_offspring, num_loci, loci_cov_matrix){
  epsilon <- MASS::mvrnorm(n_offspring, rep(0, num_loci), Sigma = loci_cov_matrix)
  select_mat <- epsilon > 0
  select_mat[] <- as.numeric(select_mat)
  matrix(rbinom(n_offspring * num_loci, 1, selection_prob) == 1,
         nrow = n_offspring,
         ncol = num_loci)
}








n_patches <- 10
dummy <- matrix(0, n_patches, n_patches)
adjacent <- abs(row(dummy) - col(dummy)) == 1
adjacent[] <- as.numeric(adjacent)


sweep(adjacent, 1, rowSums(adjacent), FUN = "/")