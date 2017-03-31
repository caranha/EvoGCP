#' Tester
#'
#' Solves a set of 3GCP problems using a specific solving method and returns the results
#'
#'
#' example of multi-core operation using the parallel library
#'
#' no_cores <- parallel::detectCores() - 1
#' cl <- parallel::makeCluster(no_cores)
#' x <- tester(pset, nfe, solver, NULL, cl = cl)
#' parallel::stopCluster(cl)
#'
#' @param problemset a list of 3GCP graphs
#' @param nfe the maximum number of evaluations per graph
#' @param solver the solver method to be used
#' @param solverpar a list with the parameters for solver
#' @return a list with summarized results
#'
#' @export
tester <- function(pset, nfe, solver, solverpar, cl = NULL, output = NULL) {

  solver_function <- paste0("solver_", tolower(solver))

  result <- if (requireNamespace("pbapply", quietly = TRUE)) {
    pbapply::pboptions(type = "timer")
    result <- pbapply::pbsapply(pset,
                                tester.solve, nfe, solver_function, solverpar,
                                cl = cl)
  } else {
    sapply(pset,
           tester.solve, nfe, solver_function, solverpar)
  }

  result <- list(violation  = unlist(result[1, ]),
                 solution   = matrix(unlist(result[2, ]), nrow = ncol(result)),
                 evaluation = unlist(result[3, ]))

  if (!is.null(output)) { save(result, file = output) }

  return (result)
}

tester.solve <- function(G, nfe, solver_function, solver_par) {
  do.call(solver_function, args = list(G, nfe, solver_par))
}

testcl <- function(pset, nfe, solver)
{
  no_cores <- parallel::detectCores() - 1
  cl <- parallel::makeCluster(no_cores)

  x <- tester(pset, nfe, solver, NULL, cl = cl)

  parallel::stopCluster(cl)
}




