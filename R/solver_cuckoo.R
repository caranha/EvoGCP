#' Cuckoo Search Solver (Toda Version)
#'
#' Solves a instance of the 3-GCP problem using the Cuckoo Search (CS)
#' implementation described in Toda et al., 2016.
#'
#' The CS algorithm begins with a random set of solutions X, and at every
#' iteration performs the following two steps:
#'
#' \itemize{
#' \item 1- For every solution x_i in X, apply the mutate.cuckoo(x_i, Levy)
#' function to generate a new solution x_u. Replace x_i with x_u if the
#' second is better.
#' \item 2- For every solution x_i in X, apply the mutate.cuckoo(x_i, Policy)
#' function with a probability \emph{pc} to generate a new solution x_u.
#' Replace x_i with x_u always (if \emph{compare} is false) or if the
#' second is better (otherwise).
#' }
#'
#' Mutate.cuckoo(x_i, policy) generates a new individual as follows:
#' an integer `m` is chosen based on the policy parameter (levy distribution,
#' uniform distribution, or fixed). Then `m`` elements from x_i are
#' changed to a random, different value.
#'
#' @param G the graph to be solved, represented by a list where G$V is the
#' number of nodes, and G$E is a |E|x2 matrix of edges.
#' @param nfe the number of function evaluations. The solver will stop after this
#' number has been exceeded.
#' @param args a list with arguments for the method. The list must contain the
#' following names:
#'
#' \itemize{
#' \item \emph{pop}: Integer > 0. The size of the solution set X
#' \item \emph{pc}: Float in [0,1]. The probability of mutation in step 2.
#' \item \emph{compare}: Boolean. Wheter the individuals in step 2 are compared
#' before being accepted.
#' \item \emph{policy}: One of "levy", "uniform", "fixed". Whether the m
#' value in step 2 is chosen from a levy distribution, uniform distribution
#' or a fixed value.
#' \item \emph{E}: Value of 'm' for \emph{policy = "fixed"}
#' \item \emph{beta}: Parameter for the levy distribution
#' \item \emph{alpha}: Parameter for the levy distribution
#' }
#'
#' @return A list with three names:
#' \itemize{
#' \item{violation}: the number of graph coloring violations of the best solution found (0 for a correct solution)
#' \item{best}: a vector with the best solution found
#' \item{evals}: the number of evaluations used by the time the solver stopped.
#' }
#'
#' @references
#' Toda Keita, Claus Aranha, Hitoshi Kanoh, "Solving the Graph Coloring Problem using Cuckoo Search",
#' Technical Report of the Information Processing Society of Japan, 2016
#'
#' @export
solver_cuckoo <- function(G, nfe, args)
{
  eval <- 0
  vio.best <- nrow(G$E) + 1
  c.best <- NULL

  # Parameters:
  pop <- args[["pop"]]                 # population size
  pc <- args[["pc"]]                   # parasitism chance (0..1)
  pcompare <- args[["compare"]]        # Perform comparison in parasitism stage
  ppolicy <- args[["policy"]]          # parasitism policy (uniform, levy, fixed)
  E <- args[["E"]]                     # fixed parasitism value
  b <- args[["beta"]]                  # beta parameter for Levy dist
  a <- args[["alpha"]]                 # alpha parameter for Levy dist

  # initial population
  P <- t(sapply(1:pop, FUN = function(x) { random_coloring(G$V) }))
  V <- apply(P, 1, evaluate, graph = G)

  vio.best <- V[order(V)[1]]
  c.best <- P[order(V)[1], ]
  eval <- pop

  while (vio.best > 0 & eval < nfe) {

    # Step 1: Perform variation on all candidates:
    newP <- cuckoo.step1(P, G, a, b)

    eval <- eval + newP$eval
    D <- (V > newP$V)           # pairwise testing children and parent
    P[D, ] <- newP$P[D, ]       # replace better children
    V[D]   <- newP$V[D]

    # Step 2: Perform variation on some candidates following "parasitism policy"
    newP <- cuckoo.step2(P, G, a, b, pc, ppolicy, E)

    eval <- eval + newP$eval
    if (pcompare == T) { D <- (newP$V < V) }        # pairwise testing children and parent
                  else { D <- (newP$V < nrow(G$E) + 1) } # If no comparison, accept all valid mutations
    P[D, ] <- newP$P[D, ]       # replace better children
    V[D]   <- newP$V[D]

    # Update best individuals
    if (min(V) < vio.best) {
      vio.best <- V[order(V)[1]]
      c.best <- P[order(V)[1], ]
    }
  }

  return(list(violation = vio.best, best = c.best, evals = eval))
}

cuckoo.mutate <- function (ind, m)
{
  v <- length(ind)
  change <- rep(0,v)
  change[sample(v,m)] <- sample(2, m, replace = T)
  ind <- ((ind - 1) + change) %% 3 + 1
  return(ind)
}

cuckoo.step1 <- function(P, G, alpha, beta) {
  M <- pmin(floor(alpha * abs(rlevy(nrow(P), beta))) + 1, G$V)
  P2 <- t(sapply(1:nrow(P), function(x) { cuckoo.mutate(P[x, ], M[x]) }))
  V2 <- apply(P2, 1, evaluate, graph = G)

  list(P = P2, V = V2, eval = nrow(P))
}

cuckoo.step2 <- function(P, G, alpha, beta, pc, ppolicy, E) {
  M <- switch(ppolicy,
              fixed = rep(E, nrow(P)),
              uniform = sample(G$V, nrow(P)),
              levy = pmin(floor(alpha * abs(rlevy(nrow(P), beta))) + 1, G$V)
  )
  M <- M * (runif(nrow(P)) <= pc)
  P2 <- t(sapply(1:nrow(P), function(x) { cuckoo.mutate(P[x, ], M[x]) }))
  V2 <- sapply(1:nrow(P), function(x) { if (M[x] > 0) { evaluate(P2[x, ], G) } else { nrow(G$E) + 1 }})

  list(P = P2, V = V2, eval = sum(M > 0))
}
