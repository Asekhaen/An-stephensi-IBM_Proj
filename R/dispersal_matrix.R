# dispersal matrix 

dispersal_frac <- 0.2    # Base dispersal rate (adjust this if necessary)

dist_matrix <- as.matrix(dist(coords, method = "euclidean"))


# create exponential dispersal kernel using the distance matrix. 
lambda <- 5 # range 

dispersal_kernel <- exp(-lambda * dist_matrix)

# set the diagonal elements to 0 to prevent self-dispersal

diag(dispersal_kernel) <- 0


# make these columns sum to 1 to get probability of moving to other patch
# *if* they left. This dispersal matrix gives the probability of the vector
# vector moving between patches
rel_dispersal_matrix <- sweep(dispersal_kernel, 1,
                              rowSums(dispersal_kernel), FUN = "/")

# sum(rel_dispersal_matrix[2,])


# normalise these to have the overall probability of dispersing to that patch,
# and add back the probability of remaining
dispersal_matrix <- dispersal_frac * rel_dispersal_matrix +
  (1 - dispersal_frac) * diag(nrow(dispersal_kernel))
