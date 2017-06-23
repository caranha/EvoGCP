#' Tester
#'
#' Solves a set of 3GCP problems using a specific solving method and returns the results
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
tester <- function(pset, nfe, solverpar, cl = NULL, output = "none") {

  solver_function <- paste0("solver_", tolower(solverpar$name))

  result <- if (requireNamespace("pbapply", quietly = TRUE)) {
    pbapply::pboptions(type = output)
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

  return (result)
}

#' Tester Batch
#'
#' Solves a large number of 3-gcp problems, with a specified set of solver and solver parameters,
#' and a specified set of graph parameters. Results are saved as the tester goes along, and it is
#' able to load previous results to avoid re-computation.
#'
#' @param solverpar.list a list where each element is a parameter list for test.solver
#' @param graphpar.list a list where each element is a parameter list for dataset
#' @param results.file a filename for saving partial and final results
#' @param graphs.file a filename for saving partial and final data sets
#' @param nfe number of maximum function evaluations for this experiment
#' @param random.seed the random number seed for this experiment
#' @param ncores the maximum number of cores to be used in this experiment
#' @return Nothing, but saves the files "results.file" and "graphs.file" with the calcultion results
#'
#' @examples
#' graphpar.list = list()
#' graphpar.list[["random90"]] = list(name = "random90", size = 5, nodes = 90, density = 3, method = "random")
#' graphpar.list[["minton90"]] = list(name = "random90", size = 5, nodes = 90, density = 3, method = "minton")
#'
#' solverpar.list = list()
#' solverpar.list[["cs.toda"]] = list(name = "cs.toda", method = "cuckoo", pop = 10, pc = 0.0001, compare = F, policy = "levy", alpha = 1, beta = 1.5, E = NULL)
#' solverpar.list[["abc.togashi"]] = list(name = "abc.togashi", method = "abc", pop = 500, onlooker = 500, c = 3, scout = 1, limit = 300)
#'
#' tester.batch(solverpar.list, graphpar.list, results.file = "foo.Rda", graphs.file = "bar.Rda", nfe = 10000, random.seed = 42, ncores = 1)
#'
#' file.remove("foo.test")
#' file.remove("bar.test")
#'
#' @export
tester.batch <- function(solverpar.list, graphpar.list, results.file, graphs.file, nfe, random.seed = 42, ncores = Inf, output = "none")
{
  assertthat::assert_that(
    is.list(solverpar.list),
    is.list(graphpar.list),
    assertthat::not_empty(c(solverpar.list, graphpar.list)),
    is.character(results.file),
    is.character(graphs.file),
    assertthat::is.count(nfe),
    is.numeric(random.seed)
  )

  graphs.list = list()
  results.list = list()

  # loading existing resources, if they exist
  if (any(file.exists(c(results.file,graphs.file)))) {
    print("graphs or result files exist. Graphs and Methods with existing names will not be re-generated.")

    if (file.exists(results.file)) load(results.file)
    if (file.exists(graphs.file)) load(graphs.file)
  }

  set.seed(random.seed)
  # loading/generating grahs
  for (g in graphpar.list)
  {
    ## Test the graphpar list
    if (!all(assertthat::has_name(g, c("name", "size", "nodes", "density", "method"))))
      stop("Malformed list of graph sets")
    if ((g["name"] %in% names(graphs.list)))
    {
      print(paste("Reading", g$name, "from file."))
    }
    else
    {
      print(paste("Generating Dataset", g$name));
      graphs.list[[g$name]] <- problemset(g$size, g$nodes, g$density, g$method)
    }
    save (graphs.list, file = graphs.file)
  }


  # Generating a random order of experiments to be executed
  set.seed(NULL)
  exp <- expand.grid(names(solverpar.list), names(graphpar.list), stringsAsFactors = F)
  exp <- exp[sample(seq(nrow(exp))), ]

  # Setup CL
  n.cores <- min(ncores, parallel::detectCores() - 1)
  cl <- if (requireNamespace("parallel", quietly = TRUE)) { parallel::makeCluster(n.cores) } else { NULL }

  for (i in seq_along(exp[, 1]))
  {
    exp.name <- paste(exp[i, 1], exp[i, 2], sep="_")
    if(exp.name %in% names(results.list))
    {
      print(paste("Reading results for method", exp[i, 1], "in dataset", exp[i, 2], "from file (",i,"out of",nrow(exp),")."))
    }
    else
    {
      print(paste("Preparing to run method", exp[i, 1], "in dataset", exp[i, 2], "(",i,"out of",nrow(exp),")."))
      set.seed(random.seed)

      par <- solverpar.list[[exp[i,1]]]
      par$name <- par$method
      results.list[[exp.name]] <- tester(graphs.list[[exp[i,2]]], nfe = nfe, solverpar = par, cl = cl, output)
    }

    save (results.list, file = results.file)
  }

  # Close CL
  if (!is.null(cl)) { parallel::stopCluster(cl) }
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
