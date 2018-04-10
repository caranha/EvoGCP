#' Firefly Solver (Fister Version)
#'
#' Solves a instance of the 3-GCP problem using the Firefly Algorithm (FFA)
#' implementation described in Fister et al., 2012.
#'
#' The FFA algorithm begins with a random set of solutions X, and at every
#' iteration performs two steps:
#'
#' \itemize{
#' \item 1- Move fireflies
#' \item 2- Heuristical Swap (optional)
#'
#' After each step, each solution is evaluated after being decoded
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
  eval <- newP$evals

  while(vio.best > 0 & eval < nfe) {
    # Step 1. Move all fireflies
    newP <- ffa.step1(W, V, G, alpha, beta, gamma, pop)
    SATUR <- newP$SATUR

    #elitism
    eval <- eval + newP$evals
    D <- (newP$V <= V)          # pairwise testing children and parent
    W[D, ] <- newP$W[D, ]       # replace better children
    P[D, ] <- newP$P[D, ]
    V[D]   <- newP$V[D]

    # Step 2. Heuristical Swap
    if(use_swap == T){
      newP <- ffa.step2(P, W, V, SATUR, G)

      #elitism
      eval <- eval + newP$evals
      D <- (newP$V <= V)          # pairwise testing children and parent
      W[D, ] <- newP$W[D, ]       # replace better children
      P[D, ] <- newP$P[D, ]
      V[D]   <- newP$V[D]
    }

    # Update best individuals
    if (min(V) <= vio.best) {
      vio.best <- V[order(V)[1]]
      c.best <- P[order(V)[1], ]
    }
  }

  return(list(violation = vio.best, best = c.best, evals = eval))

}

#move firefly
ffa.step1 <- function(W, V, G, alpha, beta, gamma, pop) {
  o <- order(V)
  W2 <- t(sapply(1:pop, ffa.move, W=W, V=V, o=o, alpha=alpha, beta=beta, gamma=gamma))

  newP <- ffa.decode(W2, G)
  return(list(P = newP$P, V = newP$V, evals = newP$eval, W=W2, SATUR=newP$SATUR))
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
  P2 <- sapply(1:nrow(W), FUN = function(x) {
    solver_dsatur(G, G$V+1, args=list(weight=W[x,], return_satur=TRUE))})

  #extrac P, V and eval
  P3 <- t(sapply(1:nrow(W), FUN = function(x) { P2[["best",x]] }))
  SATUR <- t(sapply(1:nrow(W), FUN = function(x) { P2[["satur",x]] }))
  V2 <- unlist(unname(P2["violation",]))
  eval2 <- sum(unlist(P2["evals",]))

  return(list(P = P3, V = V2, evals = eval2, SATUR = SATUR))
}

#' Heuristical Swap (step2)
#'
#' Swaps an uncolored vertex with the vertex that has highest saturation degree
#' if tie, choose randomly
#'
#' Steps:
#' 1. get the first uncolored vertex
#' 2. sort the predecessors according to the saturation degree descending
#' 3. swap the uncolored vertex with the vertex that has highest saturation degree, if tie, choose randomly (from the tie set)
#'
#' @param p vector of colors
#' @param w vector of weights
#' @param v the uncolored vertex which will be swapped
#' @param satur vector of saturation degrees
#' @return a permutation of w
ffa.heuristical_swap <- function(p, w, v, satur){
  w2 <- w

  if(!is.na(v) & v > 1){ #there is no vertex to swap, if the uncolored vertex is the first one
    o <- order(satur[1:(v-1)], decreasing = T)

    #get the vertex that has the highest satur
    U <- which(satur[1:(v-1)] == satur[o[1]]) #list of vertices with highest satur
    u <- sample(U, 1) #choose one of them randomly

    #swap v and u
    w2[v] <- w[u]
    w2[u] <- w[v]
  }

  return(w2)
}

#' Improve (step2)
#'
#' Tries to improve a solution by local search (heuristical swap)
#' Stops when improvement is detected*
#'
#' NOTE: The paper is a little bit unclear about this procedure. It says the heuristical swap operator is executed
#' until the improvements are detected, but it doens't say exactly how to proceed after a no improvement
#' Some questions of when there is no improvement:
#' should the next uncolored vertex be considered? or
#' should the next swap be on the result of the previous swap?
#'
#' This version assumes the next uncolored vertex should be considered in case of no improvement.*
#'
#'
#' @param p vector of colors
#' @param w vector of weights
#' @param v number of violations
#' @param G graph
#' @return list with the permutation of w and its related p, v and number of its functions evaluated
ffa.improve <- function(p, w, v, satur, G){
  #setup
  i <- 1
  evals <- 0
  climbing <- TRUE

  #list of uncolored vertices
  uncolored <- which(p == 0)

  #saturation degrees
  #adjacent_color <- adjacent_color_set(p, G)
  #satur <- sapply(1:ncol(adjacent_color), FUN = function(x) { sum(adjacent_color[,x]) })

  #run until improvement is detected
  while(i <= length(uncolored) & climbing){
    #move solution
    neww <- ffa.heuristical_swap(p, w, uncolored[i], satur)

    #evaluate new solution
    newp <- solver_dsatur(G, G$V+1, args=list(weight=neww, return_satur=FALSE))

    #if the new solution is better, stop. Otherwhise move again
    if(newp$violation < v){
      p <- newp$best
      w <- neww
      v <- newp$violation
      climbing <- FALSE
    }
    evals <- evals + newp$evals #+ 1
    i <- i + 1
  }

  return(list(p = p, v = v, eval = evals, w = w))
}

#Local Search
ffa.step2 <- function(P, W, V, SATUR, G){
  newP <- sapply(1:nrow(W), FUN = function(x){ ffa.improve(p=P[x,], w=W[x,], v=V[x], satur=SATUR[x,], G=G) })

  P2 <- t(sapply(1:nrow(W), FUN = function(x) { newP[["p",x]] }))
  W2 <- t(sapply(1:nrow(W), FUN = function(x) { newP[["w",x]] }))
  V2 <- unlist(newP["v",])
  eval2 <- sum(unlist(newP["eval",]))

  return(list(P = P2, V = V2, evals = eval2, W=W2))

}
