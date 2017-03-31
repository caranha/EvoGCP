
#' Random coloring generator
#'
#' Generate a random solution for the graph 3-coloring problem.
#'
#' @param N Number of vertices
#' @return A vector with N integers 1-3
random_coloring <- function(N) {
  sample(3, N, replace = T)
}

#' Connectivity
#'
#' Tests a graph for connectivity
#'
#' @param V number of nodes
#' @param E edge list represented as an |E|x2 matrix
#' @return Boolean whether the graph is connected
connectivity <- function(V,E) {

  A <- matrix(0,V,V)

  for (i in 1:V)
  {
    A[i, i] <- 1
  }

  for (i in 1:nrow(E))
  {
    A[E[i, 1], E[i, 2]] <- 1
    A[E[i, 2], E[i, 1]] <- 1
  }

  K <- A
  for (i in 1:V)
  {
    if (min(K) > 0) { return(T) }
    K <- (K > 0)
    K <- K %*% A
  }

  min(K) > 0
}

#' Random Levy
#'
#' Levy Flight distribution function, based on
#' the description by Xin-She Yang,
#' "Nature-Inspired Metaheuristic Algorithms"
#'  Second Edition, Luniver Press, (2010).
#'
#' @param N number of samples to generate
#' @param b is a scaling parameter between 0.3 and 1.99
#' @return N samples taken from the Levy Flight Distribution
rlevy <- function(N, b) {
  p.sigma <- ((gamma(1 + b) * sin(pi * b / 2)) /
               (gamma((1 + b) / 2) * b * (2 ^ ((b - 1) / 2)))) ^ (1 / b)
  q.sigma <- 1

  p <- rnorm(N, 0, p.sigma)
  q <- rnorm(N, 0, q.sigma)

  p / (abs(q) ^ (1 / b))
}

#' Similarity
#'
#' Calculate the hamming distance-based similarity between two
#' (sets of) coloring candidates.
#'
#' Led d be the number of positions with different colors between
#' coloring x1 and x2, then the similarity between x1 and x2 is
#' given by:
#'
#' 1 - d / (length of x)
#'
similarity <- function(x1, x2) {
  if (is.null(nrow(x1))) {
    d <- sum(x1 != x2)
    l <- length(x1)
  } else {
    d <- rowSums(x1 != x2)
    l <- ncol(x1)
  }

  1 - d / l
}

#' Kui's similarity
#'
#' Based on the hamming distance ''similarity'', but with added randomness
#' for extra fun. The similarity between coloring vectors x and y is
#' calculated as before, but with an extra factor R.
#'
#' 1 - R * (d / length)
#'
#' Where d is the number of positions with different colors between x and y,
#' and R is a random uniform number between 0 and 1.
#'
similarity.kui <- function(x1, x2) {
  if (is.null(nrow(x1))) {
    d <- sum(x1 != x2)
    l <- length(x1)
    r <- runif(1)
  } else {
    d <- rowSums(x1 != x2)
    l <- ncol(x1)
    r <- runif(nrow(x1))
  }

  1 - r * (d / l)
}

