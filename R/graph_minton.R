#' Generates a 3-GCP instance following the algorithm described by Minton et al.
#'
#' S. Minton, M. Johnston, A. Philips, P. Laird, "Minimizing Conflicts: a heuristic repair method
#' for constraint satisfaction and scheduling problems", Artificial Intelligence, v.58 pp.161-205, 1992
#'
#' The algorithm is defined in three steps:
#'
#' 1- Create 3 groups with N/3 nodes
#' 2- Create E edges between the groups randomly
#' 3- Accept if the graph has no unconnected components
#'
#' @param N Integer - number of vertices in the graph.
#' @param density 0..1 - edge density in the graph.
#' @return a list of pair arrays representing edges in the graph. The graph is simple and undirected.
#' @examples
#' generate_minton(10,2)
#' @export
graph_minton <- function(N, density = 2) {

  V1 <- 1:floor(N / 3)
  V2 <- (floor(N / 3) + 1):floor((N * 2) / 3)
  V3 <- (floor((N * 2) / 3) + 1):N

  en <- round(N * density)

  E12 <- expand.grid(V1, V2)
  E13 <- expand.grid(V1, V3)
  E23 <- expand.grid(V2, V3)

  Elist <- rbind(E12, E13, E23)

  if (nrow(Elist) < en) stop("Not enough vertices for a 3-colorable graph of this density")

  C <- F

  while (C == F) {
    E <- Elist[sample(nrow(Elist), en), ]
    C <- connectivity(N, E)
  }

  return(list(V = N, E = E))
}
