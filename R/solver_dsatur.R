#' DSatur Solver
#'
#' Solves the 3-GCP problem using DSatur heuristic
#'
#' DSatur heuristic function, adapted from the description by Daniel Brelaz,
#' "New methods to Color the Vertices of a Graph" Communications of the ACM. 22(4): 251–256 (1979).
#'
#' The algorithm follows 5 steps:
#' 1. Arrange the vertices by decreasing order of degrees (or any other permutation).
#' 2. Color a vertex of maximal degree with color 1.
#' 3. Choose a vertex with a maximal saturation degree.
#' 4. Color the chosen vertex with the least possible (lowest numbered) color.
#' 5. If all the vertices are colored, stop. Otherwise, return to 3.
#'
#'
#' @param G the graph to be solved, represented by a list where G$V is the number of nodes, G$E is a |E|x2 edge matrix,
#' G$adj is the adjacency list.
#' @param nfe the number of iterations allowed
#' @param args a list with arguments for the method. The list must contain the
#' following names:
#' \itemize{
#' \item \emph{weight}: List of weights, which will be used to build a permutation of vertices.
#' \item \emph{partial_solution}: List of colors. Partial solution for G.
#' \item \emph{leave_uncolored}: Boolean. If it is allowed to have uncolored vertices (then, no violation) in the solution
#' }
#'
#' @return a list containing the total number of violations of the best coloring, the best coloring
#' (a V vector of 1:3) and the total number of evaluations spent
#'
#' @examples
#' solver_dsatur(10,2)
solver_dsatur <- function(G, nfe, args){#wrapper
  assertthat::assert_that(
    all(assertthat::has_name(args,
                             c("weight", "partial_solution", "leave_uncolored")
    )))

  # Parameters:
  weight <- args[["weight"]]
  partial_solution <- args[["partial_solution"]]
  leave_uncolored <- args[["leave_uncolored"]]

  if(is.null(partial_solution)){
    partial_solution <- sample(0, G$V, replace=T)
  }

  if(leave_uncolored == T){
    result <- dsatur2(G, nfe, weight, partial_solution) #doesn't color all vertices, doesn't generate violation
  }else{
    result <- dsatur(G, nfe, weight, partial_solution) #colors all vertices, can generate violation
  }

  violation <- evaluate(result[["best"]], G)

  return(list(violation = violation, best=result[["best"]], evals=result[["evals"]]))

}

#' DSatur (colors all vertices, can generate violation)
dsatur <- function(G, nfe, weight, partial_solution){
  #Init: solution, adjacent colors of each vertex, saturation degrees
  #adjacent_color <- matrix(0, 3, G$V)
  solution <- partial_solution
  adjacent_color <- adjacent_color_set(solution)

  satur <- rep(0, G$V)
  nit <- 0

  #if tie in saturation, use this permutation
  perm <- order(weight, decreasing = T)

  while((nit < nfe) & (length(solution[solution == 0])>0)){
    #Next vertex
    v <- perm[which(satur[perm] == max(satur))[1]]

    #Lowest Numbered color
    color <- which.min(adjacent_color[,v])

    #Color vertex
    solution[v] <- color

    #update counter of adjacent colors
    adjacent_color[color, G$adj[[as.character(v)]]] <- 1

    #keep track of uncolored vertices
    uncolored_v <- which(solution==0)

    #update saturation degree (and reduce priority of already colored vertices)
    satur <- sapply(1:G$V, FUN = function(x) { if(is.element(x, uncolored_v)) {sum(adjacent_color[,x])} else -1 })

    nit <- nit + 1
  }

  return(list(best = solution, evals = nit))
}

#' DSatur (doesn't color all vertices, doesn't generate violation)
dsatur2 <- function(G, nfe, weight, partial_solution){
  #Init: solution, adjacent colors of each vertex, saturation degrees
  #adjacent_color <- matrix(0, 3, G$V)
  solution <- partial_solution
  adjacent_color <- adjacent_color_set(solution)
  satur <- rep(0, G$V)
  nit <- 0
  colored <- c() #set of colored vertices and vertices that could not be colored due to violation.

  #if tie in saturation, use this permutation
  perm <- order(weight, decreasing = T)

  while((nit < nfe) & (length(colored)<=G$V)){
    #Next vertex
    v <- perm[which(satur[perm] == max(satur))[1]]

    #Lowest Numbered color
    color <- which.min(adjacent_color[,v])

    #if it is possible to color v without violation, color it
    #otherwise, leave it with no color
    if(adjacent_color[color,v] == 0){
      #Color vertex
      solution[v] <- color

      #update counter of adjacent colors
      adjacent_color[color, G$adj[[as.character(v)]]] <- 1

      #update saturation degree
      satur <- sapply(1:G$V, FUN = function(x) { sum(adjacent_color[,x]) })
    }
    #update list of already processed vertices and reduce their priority
    colored <- c(colored, v)
    satur[colored] <- -1

    nit <- nit + 1
  }

  return(list(best = solution, evals = nit))
}

#' Auxiliar function for DSatur
#' Given a set of partial coloring of G, compute adjancent colors to each vertex
#' @param s, set of colors, size |V|
#' @return a matrix (col: vertices, row: colors) containing the colors adjacent to each vertex
#' (a V vector of 1:3) and the total number of evaluations spent
#'
#' @export
adjacent_color_set <- function(s){
  adjacent_color <- matrix(0, 3, G$V)

  for(j in 1:3){
    for(i in which(s == j)){
      adjacent_color[j, G$adj[[as.character(i)]]] <- 1
    }
  }

  return(adjacent_color)
}
