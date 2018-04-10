#' Generates a 3-GCP instance following the Culberson's Graph Generator, available at:
#' http://webdocs.cs.ualberta.ca/~joe/Coloring/Colorsrc/index.html
#'
#' Implemented methods:
#' 1. Flat graph
#' 2. Equipartite graph with independent random edge assignment
#' 3. K-Colorable graph with variability and independent random edge assignment
#'
#' Flat Graph
#' Source: http://webdocs.cs.ualberta.ca/~joe/Coloring/Generators/flat.html
#'
#' The algorithm is defined in 3 steps:
#'
#' 1- Create a permutation of nodes
#' 2- Create 3 groups with N/3 nodes
#' 3- Given two groups A and B, select edges at random from the pairs AxB, subject to:
#'  - for any a in A degree(a in AxB) <= ceil(e_{AB} / |A|) + flatness
#'  - for any b in B degree(b in AxB) <= ceil(e_{AB} / |B|) + flatness
#'  - where e_{AB} = |A|*|B|*(edge's probability)
#'
#' Discussion about this graph:
#' J. Culberson, F. Luo. “Exploring the k-Colorable landscape with iterated greedy.”
#' Cliques, Coloring, and Satisfiability DIMACS Series in Discrete Mathematics and Theoretical
#' Computer Science, Jan. 1996, pp. 245–284.
#'
#' @param N Integer - number of vertices in the graph.
#' @param probability 0..1 - edge's probability.
#' @param flatness 0.. - (flatness = 0 -> fairly uniform)
#' @param K - number of partitions (colors)
#' @return a list with two elements: V - number of nodes in the graph, and E - an edge list in the shape of a 2xE matrix. The
#' graph is guaranteed to have a valid 3-coloring
#' @examples
#' graph_culberson_flat(10,1,0)
#' @export
graph_culberson_flat <- function(N, probability=1, flatness=0) {
  #generate permutation of vertices
  Vperm <- sample(1:N)
  V1 <- Vperm[1:floor(N / 3)]
  V2 <- Vperm[(floor(N / 3) + 1):floor((N * 2) / 3)]
  V3 <- Vperm[(floor((N * 2) / 3) + 1):N]

  #generate edges
  E12 <- flat.edges(V1, V2, probability, flatness)
  E13 <- flat.edges(V1, V3, probability, flatness)
  E23 <- flat.edges(V2, V3, probability, flatness)
  Elist <- rbind(E12, E13, E23)

  #shuffle colums
  idx <- sample(1:nrow(Elist), sample(nrow(Elist), 1))
  E <- Elist
  E[idx,1] <- Elist[idx,2]
  E[idx,2] <- Elist[idx,1]
  #shuffle rows
  E <- E[sample(nrow(E)), ]

  return(list(V = N, E = E))
}

flat.edges <- function(V1, V2, probability, flatness){
  #generate random edges under the constraint:
  #for any a in A degree(a in AxB) <= ceil(e_{AB} / |A|) + flatness

  #|E|
  en <- round(length(V1)*length(V2)*probability)
  #max degree for vertices in the same group
  maxdegA <- ceiling(en/length(V1)) + flatness
  maxdegB <- ceiling(en/length(V2)) + flatness
  #to keep track of all degrees
  degree <- rep(0, length(V1)+length(V2))
  names(degree) <- c(V1, V2)

  #generate edges between groups
  E12 <- expand.grid(V1, V2)
  #shuffle
  E12 <- E12[sample(nrow(E12)),]
  #choose edges
  Eb <- sapply(1:nrow(E12), FUN = function(x) { #can it be simplified?
    e <- F
    v1 <- as.character(E12[x,][1])
    v2 <- as.character(E12[x,][2])

    if(degree[v1] <= maxdegA & degree[v2] <= maxdegB){
      degree[v1] <- degree[v1] + 1
      degree[v2] <- degree[v2] + 1
      e <- T
    }

    return(e)
  })
  E <- E12[Eb,]

  if (nrow(E) < en) stop("Cannot allocate edges for a 3-colorable graph of this edge's probability and graph's flatness")
  return(E)
}

#' Equipartite Graph (+IID)
#' Source: http://webdocs.cs.ualberta.ca/~joe/Coloring/Colorsrc/index.html
#'
#' The algorithm is defined in 3 steps:
#'
#' 1- Create a permutation of nodes
#' 1- Create 3 groups with N/3 nodes
#' 2- Assign edges by the IID method
#'
#' @param N Integer - number of vertices in the graph.
#' @param probability 0..1 - edge's probability.
#' @return a list with two elements: V - number of nodes in the graph, and E - an edge list in the shape of a 2xE matrix.  The
#' graph is guaranteed to have a valid 3-coloring, but not guaranteed to be connected.
#' @examples
#' graph_culberson_equipartite(10,1)
graph_culberson_equipartite <- function(N, probability = 1){
  if (probability < 0 | probability > 1) stop("Invalid Probability. Probability has to be in the interval [0, 1].")

  #generate permutation of vertices
  Vperm <- sample(1:N)
  V1 <- Vperm[1:floor(N / 3)]
  V2 <- Vperm[(floor(N / 3) + 1):floor((N * 2) / 3)]
  V3 <- Vperm[(floor((N * 2) / 3) + 1):N]

  #generate edges between groups
  E <- IID.edges(V1, V2, V3, probability)

  return(list(V = N, E = E))
}

#' 3-Colorable with Variability Graph (+IID)
#' Source: http://webdocs.cs.ualberta.ca/~joe/Coloring/Colorsrc/index.html
#'
#' The algorithm is defined in 3 steps:
#'
#' 1- Create a permutation of nodes
#' 2- Create 3 groups with N/3 nodes, with variability delta
#' 2- Assign edges using IID method
#'
#' @param N Integer - number of vertices in the graph.
#' @param probability 0..1 - edge's probability.
#' @param delta 0, 1, 2 - variability. It generates an uniform graph when delta = 0.
#' @return a list with two elements: V - number of nodes in the graph, and E - an edge list in the shape of a 2xE matrix. The
#' graph is guaranteed to have a valid 3-coloring, but not guaranteed to be connected.
#' @examples
#' graph_culberson_variability(10, 1, 0)
graph_culberson_variability <- function(N, probability = 1, delta = 0){
  if (probability < 0 | probability > 1) stop("Invalid Probability. Probability. has to be in the interval [0, 1].")
  if (delta < 0 | delta > 2) stop("Invalid Delta. Delta has to be in the interval [0, k-1].")

  #generate color sets
  C <- T
  while(C == T){
    if(delta > 0){
      H <- sapply(1:N, FUN = function(x) { sample(0:(delta), 1) + 1 }) #choose h at random from the range [0,delta]
      color <- sapply(H, FUN = function(x) { sample(x:3, 1) }) #for all vertices, choose the color from the range [h+1,k]
    } else{
      #delta = 0
      color <- sapply(1:N, FUN = function(x) { sample(1:3, 1) }) #for all vertices, choose the color randomly
    }

    C <- (length(unique(color)) != 3) #guarantee 3 colors
  }

  #generate permutation of vertices and assign colors
  Vperm <- sample(1:N)
  V1 <- Vperm[color == 1]
  V2 <- Vperm[color == 2]
  V3 <- Vperm[color == 3]

  #generate edges
  E <- IID.edges(V1, V2, V3, probability)

  return(list(V = N, E = E))
}

#' IID - Independent Random Edge Assignment
#' Source: https://webdocs.cs.ualberta.ca/~joe/Coloring/Generators/manual.html
#'
#' Assigns edges between nodes from different groups with a given probability
#'
#' @param V1 Vector - group of nodes with color 1.
#' @param V2 Vector - group of nodes with color 2.
#' @param V3 Vector - group of nodes with color 3.
#' @param probability 0..1 - edge's probability.
#' @return a list of edges.
IID.edges <- function(V1, V2, V3, probability){
  #generate edges between groups
  E12 <- expand.grid(V1, V2)
  E13 <- expand.grid(V1, V3)
  E23 <- expand.grid(V2, V3)

  Elist <- rbind(E12, E13, E23)

  #choose edges
  Elist <- Elist[runif(nrow(Elist)) <= probability,]

  #shuffle colums
  idx <- sample(1:nrow(Elist), sample(nrow(Elist), 1))
  E <- Elist
  E[idx,1] <- Elist[idx,2]
  E[idx,2] <- Elist[idx,1]
  #shuffle rows
  E <- E[sample(nrow(E)),]

  return(E)
}
