#' Hamming Distance PSO Solver
#'
#' Solves a instance of the 3-GCP problem using the Hamming Distance PSO (HDPSO)
#' implementation described in Aoki et al., 2015.
#'
#' The HDPSO algorithm begins with a random set of solutions P. It keeps a set of
#' best solution found at each position of X (P_best), a set of solutions
#' found at the previous iteration (P_prev), and the best solution found so far
#' (g_best).
#'
#' The algorithm define the similarity between two solutions, x and y, as:
#' \code{similarity(x,y) = 1 - (sum(x != y) / length(x))}
#'
#' At each iteration, it generates a new set of solution based on the current
#' set of solutions. For each solution x_i in P, it performs the following steps:
#' \itemize{
#' \item 1- Calculate Vrand = w * similarity(x_i, P_prev_i), Vpbest =
#' c1 * r1 * similarity(x_i, P_best_i), Vgbest = c2 * r2 * similarity(x_i, g_best)
#' \item 2- w, c1, and c2 are parameters, r1 and r2 are random numbers taken
#' from [0,1] every iteration
#' \item 3- For each k element of x_i, selects a new element using one of three options:
#' a random color, the kth color from P_best_i, the kth color from g_best.
#' \item 4- The choice at the previous step is random, but weighted by Vrand,
#' Vpbest and Vgbest, respectively.
#' }
#'
#' After the new set of solutions is generated, the current set replaces
#' the previous set, and the new set replaces the current set. The new set
#' is evaluated, and P_best and G_best are updated as necessary.
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
#' \item \emph{w}: Float in [0,1]. Weight for randomly choosing new colors.
#' \item \emph{c1}: Float in [0,1]. Weight for choosing a color from Pbest.
#' \item \emph{c2}: Float in [0,1]. Weight for choosing a color from Gbest.
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
#' Takuya Aoki, Claus Aranha, Hitoshi Kanoh, "PSO Algorithm with Transition Probability Based on Hamming Distance
#' for Graph Coloring Problem", IEEE International Conference On Systems, Man and Cybernetics, 2015.
#'
#' @export
solver_hdpso <- function(G, nfe, args)
{
  eval <- 0
  vio.best <- nrow(G$E) + 1
  c.best <- NULL

  testthat::expect_named(args,
                         c("pop", "w", "c1", "c2"),
                         ignore.order = TRUE)

  # Parameters:
  pop <- args[["pop"]]                 # population size
  w <- args[["w"]]
  c1 <- args[["c1"]]
  c2 <- args[["c2"]]

  # initial random population (P will also be the local best)
  P <- t(sapply(1:pop, FUN = function(x) { random_coloring(G$V) }))
  P_best <- P
  P_prev <- P

  V_best <- apply(P_best, 1, evaluate, graph = G)

  vio.best <- V_best[order(V_best)[1]]
  c.best <- P_best[order(V_best)[1], ]
  eval <- pop

  while (vio.best > 0 & eval < nfe) {
    # Step 1. Apply mutation to all individuals. Replace if better, increase age if not.
    newP <- hdpso.update(P, G, P_prev, P_best, c.best, w, c1, c2)
    eval <- eval + newP$eval

    # updating previous population
    P_prev <- P
    P <- newP$P

    # updating pbest and gbest
    D     <- (newP$V <= V_best)
    P_best[D,] <- newP$P[D, ]
    V_best[D]  <- newP$V[D]

    # Updating Gbest (testing if new best solution is better than current best)
    if (min(V_best) <= vio.best) {
      vio.best <- V_best[order(V_best)[1]]
      c.best <- P_best[order(V_best)[1], ]
    }
  }

  return(list(violation = vio.best, best = c.best, evals = eval))
}

hdpso.update <- function(P, G, P_prev, P_best, g_best, w, c1, c2)
{
  r <- runif(2)
  G_best <- matrix(g_best, nrow(P), ncol(P), byrow = T)

  # Choice probabilities
  Vrand = w * similarity(P, P_prev)
  Vpbest = c1 * r[1] * similarity(P, P_best)
  Vgbest = c2 * r[2] * similarity(P, G_best)

  Vrand[Vrand + Vpbest + Vgbest == 0] <- 1 # if sum of probs is 0, make it random.

  # New individuals are composed of random colors, colors from P_best or G_best, depending
  # on the probabilities Vrand, Vpbest, Vgbest

  P2 <- hdpso.mutate(Vrand, Vpbest, Vgbest, P_best, G_best)
  V2 <- apply(P2, 1, evaluate, graph = G)
  list(P = P2, V = V2, eval = nrow(P2))
}

hdpso.mutate <- function(Vrand, Vpbest, Vgbest, P_best, G_best)
{
  P_rand <- matrix( sample(3, prod(dim(P_best)), replace = T),
                    nrow(P_best), ncol(P_best))
  sel <- t(sapply(seq_along(Vrand),
            function(x) { sample(3, ncol(P_best), replace = T, prob = c(Vrand[x], Vpbest[x], Vgbest[x])) }
          ))
  P_rand[sel == 2] <- P_best[sel == 2]
  P_rand[sel == 3] <- G_best[sel == 3]
  P_rand
}


