
#' Randomly generates a 3-GCP instance
#'
#' @param N Integer - number of vertices in the graph.
#' @param density 0..1 - edge density in the graph.
#' @return a list of pair arrays representing edges in the graph. The graph is simple and undirected.
#' @examples
#' generate_random(10,2)
graph_random <- function(N, density = 2) {

  V <- 1:N
  en <- round(N * density)

  E <- t(sapply(1:en, function(x) sort(sample(V, 2))))
  E <- E[order(E[, 1], E[, 2]), ]
  E <- E[!duplicated(E), ]

  while (nrow(E) < en) {
    Ea <- t(sapply(1:(en - nrow(E)), function(x) sort(sample(V, 2))))
    E <- rbind(E, Ea)

    E <- E[order(E[, 1], E[, 2]), ]
    E <- E[!duplicated(E), ]
  }

  return(list(V = N, E = E))
}

