#' Generates a random, valid 3-GCP instance, without connectivity guarantee
#'
#' @param N Integer - number of vertices in the graph.
#' @param density 0..1 - edge density in the graph.
#' @return a list with two elements: V - number of nodes in the graph, and E - an edge list in the shape of a 2xE matrix. The
#' graph is guaranteed to have a valid 3-coloring, but not guaranteed to be connected.
#' @examples
#' generate_random(10,2)
graph_random <- function(N, density = 2) {

  V1 <- 1:floor(N / 3)
  V2 <- (floor(N / 3) + 1):floor((N * 2) / 3)
  V3 <- (floor((N * 2) / 3) + 1):N

  en <- round(N * density)

  E12 <- expand.grid(V1, V2)
  E13 <- expand.grid(V1, V3)
  E23 <- expand.grid(V2, V3)

  Elist <- rbind(E12, E13, E23)
  if (nrow(Elist) < en) stop("Not enough vertices for a 3-colorable graph of this density")

  E <- Elist[sample(nrow(Elist), en), ]

  return(list(V = N, E = E))
}

