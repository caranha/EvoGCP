#' Artificial Bee Colony Solver (Togashi Version)
#'
#' Solves a instance of the 3-GCP problem using the Artificial Bee Colony (ABC)
#' implementation described in Togashi et al., 2017.
#'
#' The ABC algorithm begins with a random set of solutions X, and at every
#' iteration performs the following three steps:
#'
#' \itemize{
#' \item 1- For every solution x_i in X, find x_j (\eqn{i != j}), and
#' calculates x_u using `mutate.abc(x_i, x_j)`. Evaluate x_u and replace x_i if better.
#' \item 2- Select \eqn{n = onlooker} solutions x_i from X (with repetition), with
#' probability proportional to their fitness. Apply Step 1 on these solutions.
#' \item 3- Select \eqn{n = scout} solutions x_i from X where AGE(x_i) > \eqn{limit}.
#' replace x_i with a random solution.}
#'
#'
#' Mutate.abc(x_i, x_j, c) generates a new individual as follows: `c` elements
#' are choosen randomly from x_j, and copied into x_i.
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
#' \item \emph{onlooker}: Integer > 0. The number of solutions chosen in step 2 of the algorithm.
#' \item \emph{scout}: Integer > 0. The number of solutions chosen in step 3 of the algorithm.
#' \item \emph{limit}: Integer > 0. Minimum number of iterations without improvement before a
#' solution will be considered for step 3.
#' \item \emph{c}: Number of elements in a solution exchanged during a mutation step.
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
#' Yuuya Togashi, Claus Aranha, Hitoshi Kanoh, "Artificial Bee Colony Algorithm for Solving
#' Graph Coloring Problem", Proceedings of the IPSJ, 2017
#'
#' @export
solver_abc <- function(G, nfe, args)
{
  eval <- 0
  vio.best <- nrow(G$E) + 1
  c.best <- NULL

  assertthat::assert_that(
    all(assertthat::has_name(args,
                             c("pop", "onlooker", "scout", "limit", "c")
    )))


  # Parameters:
  pop <- args[["pop"]]                 # population size
  pop.onlooker <- args[["onlooker"]]   # repetitions for step 2
  pop.scout <- args[["scout"]]         # repetitions for step 3
  limit <- args[["limit"]]             # minimum age value for step 3 eligibility
  c <- args[["c"]]                     # parameter for the mutation operator

  # initial population
  P <- t(sapply(1:pop, FUN = function(x) { random_coloring(G$V) }))
  V <- apply(P, 1, evaluate, graph = G)
  A <- rep(0, pop)

  vio.best <- V[order(V)[1]]
  c.best <- P[order(V)[1], ]
  eval <- pop

  while (vio.best > 0 & eval < nfe) {
    # Step 1. Apply mutation to all individuals. Replace if better, increase age if not.
    newP <- abc.step1(P, G, c)

    eval <- eval + newP$eval
    D     <- (newP$V <= V)   # pairwise testing children and parent
    P[D,] <- newP$P[D, ]    # replace better children
    V[D]  <- newP$V[D]
    A[D]  <- 0              # reset new solutions
    A[!D] <- A[!D]+1        # increase age of old solutions

    # Step 2. Apply mutation to #onlookerbee individuals, choosen by 1-vio/viomax proportion. Replace if better
    #         Note, choice probability does not update here. Age does not update here.
    newP <- abc.step2(P, G, V, c, pop.onlooker)

    eval <- eval + newP$eval
    D      <- (newP$V <= V)             # replace all better individuals
    P[D, ] <- newP$P[D, ]              # replace better children
    V[D]   <- newP$V[D]
    A[D]   <- 0

    # Step 3. For all individuals with age > limit, sample up to #scoutbee individuals
    newP <- abc.step3(P, G, A, limit, pop.scout)

    eval <- eval + newP$eval
    D <- (newP$V < nrow(G$E) + 1)      # Replace all valid individuals
    P[D, ] <- newP$P[D, ]              # replace better children
    V[D]   <- newP$V[D]
    A[D]   <- 0

    # Update best individuals
    if (min(V) <= vio.best) {
      vio.best <- V[order(V)[1]]
      c.best <- P[order(V)[1], ]
    }
  }

  return(list(violation = vio.best, best = c.best, evals = eval))
}

abc.step1 <- function(P, G, c) {
  M <- sapply(1:nrow(P), function(x) { sample((1:nrow(P))[-x], 1) }) # selects k for every individual
  P2 <- t(sapply(1:nrow(P), function(x) { abc.mutate(P[x,], P[M[x], ], c) }))
  V2 <- apply(P2, 1, evaluate, graph = G)

  list(P = P2, V = V2, eval = nrow(P2))
}

abc.step2 <- function(P, G, V, c, pop.onlooker) {

  sel.prob <- 1 - (V / (nrow(G$E) + 1))
  M1 <- sample(1:nrow(P), pop.onlooker, replace = T, prob = sel.prob)
  M2 <- sapply(M1, function(x) { sample((1:nrow(P))[-x], 1) })
  P2 <- t(sapply(1:length(M1), function(x) { abc.mutate(P[M1[x],], P[M2[x], ], c) }))
  V2 <- apply(P2, 1, evaluate, graph = G)

  o <- order(V2)
  idx <- o[!duplicated(M1[o])]
  Mi <- M1[idx]
  Vi <- V2[idx]
  Pi <- P2[idx, ]

  P[Mi, ] <- Pi
  V[Mi] <- Vi

  list(P = P, V = V, eval = pop.onlooker)
}

abc.step3 <- function(P, G, A, limit, pop.scout) {

  M <- seq(nrow(P))[(A > limit)]
  M <- head(M[order(A[M], decreasing = T)],pop.scout)
  V <- array(nrow(G$E) + 1, nrow(P))
  P[M, ] <- t(vapply(seq_along(M), FUN = function(x) { random_coloring(G$V) }, seq.int(G$V)))
  V[M]   <- apply(P[M, , drop = F], 1, evaluate, graph = G)
  list(P = P, V = V, eval = length(M))
}

abc.mutate <- function(ind, ind.k, c)
{
  m <- sample(length(ind), c)
  ind[m] <- ind.k[m]

  return(ind)
}
