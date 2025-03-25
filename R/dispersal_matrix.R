# build a function for distance matrix based of on the population x & y coords.

create_dist_matrix <- function(pop) {
  # This ensures we get the coordinates of the patches
  coords <- pop |>
    bind_rows() |>
    select(patch, x, y) |>
    distinct()            
  # Calculate pairwise Euclidean distances between patches
  dist_matrix <- as.matrix(dist(coords[, c("x", "y")], method = "euclidean"))
  # Return the distance matrix
  return(dist_matrix)
}


#create distance matrix using function
dist_matrix <- create_dist_matrix(grown_pop)

# create exponential dispersal kernel using the distance matrix. 
lambda <- 1/mean(dist_matrix) # range is = 1/average distance in the matrix

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