#' Random Solver
#'
#' Solves the 3-GCP problem by generating random colorings
#'
#' @param G the graph to be solved, represented by a list where G$V is the number of nodes, and G$E is a |E|x2 edge matrix.
#' @param nfe the number of function evaluations allowed
#' @param args a list of arguments -- none is used in this particular method
#' @return a list containing the total number of violations of the best coloring, the best coloring
#' (a V vector of 1:3) and the total number of evaluations spent
#'
#' @export
solver_random <- function(G, nfe, args)
{
  eval <- 0
  vio <- Inf

  vio.best <- Inf
  c.best <- NULL

  while (vio > 0 & eval < nfe) {
    c <- random_coloring(G$V)
    vio <- evaluate(c,G)
    eval <- eval + 1
    if (vio < vio.best) {
      vio.best <- vio
      c.best <- c
    }
  }

  return(list(violation = vio.best, best = c.best, evals = eval))
}
