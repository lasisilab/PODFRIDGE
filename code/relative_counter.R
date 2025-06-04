count_first_degree <- function(n_sibs, n_child) {
  # Calculate the number of first-degree relatives
  # First-degree relatives include parents, siblings, and children
  n_parents <- 2  # Assuming two parents
  n_first_degree <- n_sibs + n_child + n_parents
  return(n_first_degree)
}

count_second_degree <- function(n_half_sibs, n_nephews_nieces, n_aunts_uncles) {
  # Calculate the number of second-degree relatives
  # Second-degree relatives include aunts, uncles, nephews, nieces and half-siblings
  n_second_degree <- n_half_sibs + n_aunts_uncles + n_nephews_nieces
  return(n_second_degree)
}

count_third_degree <- function(n_cousins) {
  # Calculate the number of third-degree relatives
  # Third-degree relatives include cousins
  n_third_degree <- n_cousins
  return(n_third_degree)
}
