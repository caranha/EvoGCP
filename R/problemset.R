
#' Problemset Generator
#'
#' Generate a set of 3GCP problems
#'
#' @param P number of problems to be generated
#' @param N number of vertices in each problem
#' @param density edge density in each problem.
#' @param method a string defining the function for problem generation. Will load a "graph_method" from the package
#' @param x addictional parameter for some specific methods (e.g. flatness, variability etc)
#' @return a list with P instances of 3GCP problems with N vertices and N*density edges
#'
#' @export
problemset <- function(P, N, density, method, x=NULL)
{
  function_name <- paste0("graph_", tolower(method))

  if(is.null(x)){
    lapply(1:P, function(t) do.call(function_name, args = list(N, density)))
  } else{
    lapply(1:P, function(t) do.call(function_name, args = list(N, density, x)))
  }
}
