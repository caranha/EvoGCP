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
#' \item \emph{lb}: Integer > 0. The lower bound of the weight set
#' \item \emph{ub}: Integer > 0. The upper bound of the weight set
#' \item \emph{alpha}: Integer > 0. The size of the randomness move withim the search space
#' \item \emph{beta}: Integer > 0. The attractiveness at distance 0
#' \item \emph{gamma}: Integer > 0. The absorbtion coefficient
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
solver_ffa <- function(G, nfe, args)
{
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

  # Initial population
  W <- t(sapply(1:pop, FUN = function(x) { random_weight(G$V, lb, ub) }))
  P <- DSatur(W, G)
  V <- apply(P, 1, evaluate, graph = G)

  vio.best <- V[order(V)[1]]
  c.best <- P[order(V)[1], ]
  eval <- pop

  while (vio.best > 0 & eval < nfe) {
    # Step 1. Move all fireflies
    newW <- ffa.step1(W, P, V, G, alpha, beta, gamma)

    eval <- eval + newW$eval
    #elitism() todo
    W <- newW$W
    P <- newW$P
    V <- newW$V

    # Update best individuals
    if (min(V) <= vio.best) {
      vio.best <- V[order(V)[1]]
      c.best <- P[order(V)[1], ]
    }
  }

  return(list(violation = vio.best, best = c.best, evals = eval))

}

ffa.step1 <- function(W, P, V, G, alpha, beta, gamma) {
  o <- order(V)
  W2 <- t(sapply(1:pop, ffa.move, W=W, V=V, o=o, alpha=alpha, beta=beta, gamma=gamma))
  P2 <- DSatur(W2, G)
  V2 <- apply(P2, 1, evaluate, graph = G)

  list(W = W2, P = P2, V = V2, eval = nrow(P2))
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

ffa.move2 <- function(i){##not done
  ##same as step1+move, but more compact
  #w <- current position + sum(attractions) + random step

  M <- sapply(1:pop, FUN = function(x) { seq(nrow(W))[(V > V[x])] }) #todo: optimize
  M_ind <- seq(pop)[sapply(1:pop, FUN = function(x) { length(M[[x]])>0 })]

  #r = distance between all attracted fireflies to its attracting ones
  R <- sapply(M_ind, FUN = function(x) { sapply(M[[x]], FUN = function(y) { Euclidean.distance(W[x,], W[y,]) }) })
  #attraction = beta*exp(-gamma*r^2)(wi-wj)
  A1 <- sapply(1:length(M_ind), FUN = function(x) { sapply(R[[x]], FUN = function(r) { beta*exp(-gamma*(r^2)) }) })
  A2 <- sapply(M_ind, FUN = function(x) { t(sapply(M[[x]], FUN = function(y) { W[x,]-W[y,] })) })
  #sum of attractions
  A <- sapply(1:length(M_ind), FUN = function(x) { A1[[x]]*A2[[x]] })
  A <- t(sapply(1:length(M_ind), FUN = function(x) { sapply(1:G$V, FUN = function(y) { sum(A[[x]][,y]) }) }))
  #random step of each firefly
  RND <- sapply(1:length(M_ind), FUN = function(x) { alpha*(runif(1)-(1/2)) })
  #w <- current position + sum(attractions) + random step
  W2 <- t(sapply(1:length(M_ind), FUN = function(x) { W[M_ind[x],] + A[x,] + RND[x] }))

  #todo: build W

  return(w)
}

#ffa.step2() <- function(){
ffa.heuristical_swap <- function(){

}

#**INCOMPLETE todos:
#1. how use the weights?
#2. once this function is done, compute the adjacency list outside this function
#3. what does <least possible (lowest numbered) color> mean?  less used color in the graph OR less used color for that vertex OR ??
#4. what if all 3 colors have been used already? should leave the vertex with no color? can this case actually happen?
DSatur <- function(w, G){
  #D. Brelaz: New methods to color vertices of a graph. Communications of the ACM. 22(4): 251–256 (1979).

  p <- rep(0, G$V)
  adj <- tapply(c(G$E[[2]],G$E[[1]]), c(G$E[[1]],G$E[[2]]), unique)
  c <- matrix(0, 3, G$V) #adjacent colors of each vertex

  #1. Arrange the vertices by decreasing order of degrees.
  d <- table(unlist(G$E))
  o <- order(d, decreasing = T)

  #2. Color a vertex of maximal degree with color 1. (todo: put all lines below inside the loop)
  v <- as.numeric(names(d[o[1]])) #vertex to color
  color <- 1

  p[v] <- 1
  c[color, adj[[names(d[o[1]])]]] <- c[color, adj[[names(d[o[1]])]]] + 1 #update number of times each color has been used
  rho <- sapply(1:G$V, FUN = function(x) { sum(c[,x]) }) #update saturation degrees of the uncolored vertices

  #5. If all the vertices are colored, stop. Otherwise, return to 3.
  while(length(p[p == 0])>0){
    #3. Choose a vertex with a maximal saturation degree.
    v <- which.max(rho)

    #4. Color the chosen vertex with the least possible (lowest numbered) color.
    color <- which.min(c[,v])
    p[v] <- color #color
    c[color, adj[[v]]] <- c[color, adj[[v]]] + 1
    uncolored_v <- which(p==0, arr.ind = TRUE) #vertices with no color
    rho <- sapply(1:G$V, FUN = function(x) { if(is.element(x, uncolored_v)) {sum(c[,x])} else -1 }) #update saturation degree
  }
  return(p)

}

euclidean_distance <- function(a, b){
  if (length(a) == length(b))
    return(sqrt(sum(sapply(1:length(a), FUN = function(x) { (a[x] - b[x])^2 }))))
  return(NULL)
}

random_weight <- function(N, lb, ub){
  return(sapply(1:N, FUN = function(x) { (ub - lb) *  runif(1) + lb }))
}
