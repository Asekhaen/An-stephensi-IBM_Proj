
# random positions of target loci
n_loci_total <- 100
n_target_loci <- 20
target_loci_position <- sort(sample.int(n_loci_total,
                                        n_target_loci,
                                        replace = FALSE))

target_loci_position

# exponential model for correlated random effects

# range controls how the correlation between two loci (in logit probability)
# depends on how many positions separate two loci
range <- 50
# var controls how variable different loci are (in logit probability)
var <- 0.5

# define the covariance function on logit probabilities across the target loci
distmat <- as.matrix(dist(target_loci_position)) ^ 2
covar <- var * exp(-distmat / range)

# set the mean probability of a thing hapening at each locus
mean_prob <- 0.5

# simulate correlated probabilities of things happening
epsilon <- MASS::mvrnorm(1, rep(0, n_target_loci), Sigma = covar)
prob <- plogis(qlogis(mean_prob) + epsilon)
selection <- rbinom(n_target_loci, 1, prob)
plot(selection ~ target_loci_position,
     ylim = c(0, 1))
abline(h = mean_prob,
       lty = 2)

