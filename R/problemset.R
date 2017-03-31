
#' Problemset Generator
#'
#' Generate a set of 3GCP problems
#'
#' @param P number of problems to be generated
#' @param N number of vertices in each problem
#' @param density edge density in each problem
#' @param method a string defining the function for problem generation. Will load a "graph_method" from the package
#' @return a list with P instances of 3GCP problems with N vertices and N*density edges
#'
#' @export
problemset <- function(P, N, density, method)
{

  function_name <- paste0("graph_", tolower(method))
  lapply(1:P, function(t) do.call(function_name, args = list(N, density)))

}
