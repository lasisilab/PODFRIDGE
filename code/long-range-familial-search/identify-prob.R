
# number of children simulator based on zero-inflated negative binomial distribution
number_of_children_zinb <- function(mu, size, inflation) {
  # Generate a random number of children based on zero-inflated negative binomial distribution
  if (runif(1) < inflation) {
    return(0)  # With probability 'inflation', return 0 children
  } else {
    return(rnbinom(1, size = size, mu = mu))  # Otherwise, return a random number of children
  }
}

# calculate the number of the (g-1)-th cousins (view siblings as the 0-th cousins) with family size r
count_cousins <- function(g, r) {
  k_ancestor_couple = 2^(g-1)
  k_cousins = k_ancestor_couple*(r[g]-1) # remove the descendants of the p-th ancestor
  if(g>1) {
    for(i in 1:(g-1)) {
      k_cousins = k_cousins * r[i]
    }
  }
  return(k_cousins)
}


# family size simulator
family_size_simulator <- function(n, type = "ZINB", mu = 3, size = 1, inflation = 0.06) {
  family_size <- array(NA, c(n, 5))
  colnames(family_size) <- c("number of half-siblings",
                             "number of siblings",
                             "number of first degree cousins",
                             "number of second degree cousins",
                             "number of third degree cousins")
  family_size[, 1] <- 0 # number of half-siblings are set to 0 in this simulation

  gmax = 4 # maximum number of generations to consider
  for (i in 1:n) {
    num_child_vector <- numeric(4)
    for (g in 1:gmax) {
      if(type == "ZINB") {
        num_child_vector[g] <- number_of_children_zinb(mu, size, inflation)
      } else {
        stop("Unknown type of family size simulator")
      }
    }

    for (g in 1:gmax) {
      r = num_child_vector[1:g]
      # Ensure that the number of children of the g-th generation is not zero
      while(r[g] == 0) {
        if(type == "ZINB") {
          r[g] <- number_of_children_zinb(mu, size, inflation)
        } else {
          stop("Unknown type of family size simulator")
        }
      }
      family_size[i, g + 1] <- count_cousins(g, r)
    }

  }

}
