
# number of children simulator based on zero-inflated negative binomial distribution
number_of_children_zinb <- function(mu, size, inflation) {
  # Generate a random number of children based on zero-inflated negative binomial distribution
  if (runif(1) < inflation) {
    return(0)  # With probability 'inflation', return 0 children
  } else {
    return(rnbinom(1, size = size, mu = mu))  # Otherwise, return a random number of children
  }
}

# calculate the probability of observing a match
genetic_match <- function(g, m = 6, num_chrs = 22, genome_size = 35, min_num_seg = 2) {
  m = m / 100 # Convert m from cM to fraction
  f = exp(-2 * g * m) / 2^(2 * g - 2) # calculate the probability of sharing a detectable IBD segment
  pr = 1 - pbinom(min_num_seg - 1, num_chrs + genome_size * 2 * g, f) # calculate probability
  return(pr)
}

# family size simulator
family_size_simulator <- function(n, type = "ZINB", mu = 3, size = 3, inflation = 0.1) {
  family_size <- array(NA, c(n, 5))
  colnames(family_size) <- c("n_sibs",
                             "n_half_sibs",
                             "n_first_degree",
                             "n_second_degree",
                             "n_third_degree")
  family_size[, 2] <- 0 # number of half-siblings are set to 0 in this simulation

  for(i in 1:n) {
    n_parents <- 2 # Assuming two parents
    if(type == "ZINB") {
      n_child <- number_of_children_zinb(mu, size, inflation)  # Generate number of children using ZINB
      n_sibs <- -1
      while(n_sibs < 0) {
        n_sibs <- number_of_children_zinb(mu, size, inflation) - 1  # Generate number of siblings using ZINB
      }
      n_first_degree <- n_sibs + n_child + n_parents  # First-degree relatives include parents, siblings, and children
      family_size[i, 1] <- n_sibs
      family_size[i, 3] <- n_first_degree

      n_half_sibs <- 0  # Set number of half-siblings to 0
      n_aunts_uncles <- 0
      for(j in 1:2) {
        n_tem <- -1
        while(n_tem < 0) {
          n_tem <- number_of_children_zinb(mu, size, inflation) - 1  # Generate number of aunts/uncles using ZINB
        }
        n_aunts_uncles <- n_aunts_uncles + n_tem  # Sum the number of aunts/uncles
      }
      n_nephews_nieces <- 0
      if(n_sibs > 0) {
        for(j in 1:n_sibs) {
          n_nephews_nieces <- n_nephews_nieces + number_of_children_zinb(mu, size, inflation)  # Generate number of nephews/nieces using ZINB
        }
      }
      n_second_degree <- n_half_sibs + n_aunts_uncles + n_nephews_nieces  # Second-degree relatives include aunts, uncles, nephews, nieces and half-siblings
      family_size[i, 4] <- n_second_degree

      n_cousins <- 0
      if(n_aunts_uncles > 0) {
        for(j in 1:n_aunts_uncles) {
          n_cousins <- n_cousins + number_of_children_zinb(mu, size, inflation)  # Generate number of cousins using ZINB
        }
      }
      family_size[i, 5] <- n_cousins
    } else {
      stop("Invalid type. Use 'ZINB' or 'Poisson'.")
    }
  }

  family_size_df = data.frame(family_size)
  return(family_size_df)
}
