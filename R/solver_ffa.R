#' Firefly Solver (Fister Version)
#'
#' Solves a instance of the 3-GCP problem using the Firefly Algorithm (FFA)
#' implementation described in Fister et al., 2012.
#'
#' The FFA algorithm begins with a random set of solutions X, and at every
#' iteration performs the following three*** steps:
#'
#' \itemize{
#' \item 1- ***
#' \item 2- ***
#' \item 3- ***
#'
#'
#' Move.ffa() ***
#'
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
#' \item \emph{lb}: Number, 0 <= lv <= ub. The lower bound of the values in the weight set
#' \item \emph{ub}: Number, 0 <= lv <= ub. The upper bound of the values in the weight set
#' \item \emph{alpha}: Integer > 0. The size of the randomness move withim the search space
#' \item \emph{beta}: Integer > 0. The attractiveness at distance 0
#' \item \emph{gamma}: Integer > 0. The absorbtion coefficient
#' \item \emph{swap}: Boolean. If it is to use the heuristical swap.
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
#' I. Fister Jr., X.-S. Yang, I. Fister, J. Brest, Memetic firefly algorithm for combinatorial
#' optimization, in Bioinspired Optimization Methods and their Applications (BIOMA 2012),
#' B. Filipic and J.Silc, Eds. Jozef Stefan Institute, Ljubljana, Slovenia, 2012
#'
#' @export
solver_ffa <- function(G, nfe, args){
  eval <- 0
  vio.best <- nrow(G$E) + 1
  c.best <- NULL

  assertthat::assert_that(
    all(assertthat::has_name(args,
                             c("pop", "lb", "ub", "alpha", "beta", "gamma")
    )))

  # Parameters:
  pop <- args[["pop"]]              # population size
  lb <- args[["lb"]]                # weight's lower bound
  ub <- args[["ub"]]                # weight's upper bound
  alpha <- args[["alpha"]]          # size of the randomness move withim the search space
  beta <- args[["beta"]]            # attractiveness at distance 0
  gamma <- args[["gamma"]]          # absorbtion coefficient
  use_swap <- args[["swap"]]        # heuristical swap?


  #DSatur heuristic setup
  G$adj <- adjacency_list(G) #build adjacency list

  # Initial population
  W <- t(sapply(1:pop, FUN = function(x) { random_weight(G$V, lb, ub) }))

  #decode
  newP <- ffa.decode(W, G)
  P <- newP$P
  V <- newP$V

  vio.best <- V[order(V)[1]]
  c.best <- P[order(V)[1], ]
  eval <- newP$eval

  while(vio.best > 0 & eval < nfe) {
    # Step 1. Move all fireflies
    newW <- ffa.step1(W, V, alpha, beta, gamma, pop)

    # Decode and evaluate
    newP <- ffa.decode(W, G)

    #elitism
    eval <- eval + newP$eval
    D <- (newP$V <= V)          # pairwise testing children and parent
    W[D, ] <- newW[D, ]       # replace better children
    P[D, ] <- newP$P[D, ]
    V[D]   <- newP$V[D]

    # Update best individuals
    if (min(V) <= vio.best) {
      vio.best <- V[order(V)[1]]
      c.best <- P[order(V)[1], ]
    }

    # Step 2. Heuristical Swap (not implemented yet)
    use_swap <- F
    if(use_swap == T){
      newW <- ffa.heuristical_swap(NULL,NULL,NULL)

      # Decode and evaluate
      newP <- ffa.decode(W, G)

      # ... #
    }
  }

  return(list(violation = vio.best, best = c.best, evals = eval))

}

#move firefly
ffa.step1 <- function(W, V, alpha, beta, gamma, pop) {
  o <- order(V)
  W2 <- t(sapply(1:pop, ffa.move, W=W, V=V, o=o, alpha=alpha, beta=beta, gamma=gamma))
  return(W2)
}

ffa.move <- function(i, W, V, o, alpha, beta, gamma){
  w <- W[o[i],]

  M <- seq(i, nrow(W))[(V[o[i:nrow(W)]] > V[o[i]])]
  if(length(M) > 0){
    #distance
    R <- sapply(M, FUN = function(x) { euclidean_distance(W[o[i],], W[o[x],]) })
    #attraction
    A1 <- sapply(1:length(M), FUN = function(x) { beta*exp(-gamma*(R[x]^2)) })
    A2 <- t(sapply(M, FUN = function(x) { W[o[x],]-W[o[i],] }))
    A <- sapply(1:ncol(W), FUN = function(x) { sum((A1*A2)[,x]) })
    #random step
    rnd <- alpha*(runif(1)-(1/2))

    #w <- current position + sum(attractions) + random step
    w <- W[o[i],] + A + rnd
  }
  return(w)
}

ffa.decode <- function(W, G){
  #decode
  P2 <- sapply(1:pop, FUN = function(x) {
    solver_dsatur(G, G$V+1, args=list(weight=W[x,], partial_solution=NULL, leave_uncolored=F))})

  #extrac P, V and eval
  P3 <- t(sapply(1:pop, FUN = function(x) { P2[["best",x]] }))
  V2 <- unlist(P2["violation",])
  eval2 <- sum(unlist(P2["evals",]))

  return(list(P = P3, V = V2, eval = eval2))
}

#' Heuristical Swap
#'
#' Swaps an uncolored vertex with the vertex that has highest saturation degree
#' if tie, choose randomly
#'
#' @param p vector of colors
#' @param w vector of weights
#' @param rho vector of saturation degrees
#' @return a permutation of w
ffa.heuristical_swap <- function(p, w, rho){
  #runs till improvement?
  #1. get the first uncolored vertex
  #2. sort the predecessors according to the saturation degree descending
  #3. swap the uncolored vertex with the vertex that has highest saturation degree, if tie, choose randomly (from the tie set)

  # SKETCH (not implemented) #
  w2 <- w

  #Uncolored vertex
  v1 <- which(p == 0)[1]

  #2. todo
  #rho <- sapply(1:length(G), FUN = function(x) { sum(c[,x]) })
  #o <- order()
  v2 <- 1

  #3.
  w2[v1] <- w[v2]
  w2[v2] <- w[v1]

  return(w2)

}
