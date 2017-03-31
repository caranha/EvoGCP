
evaluate <- function(colors, graph)
{
  edges = graph$E
  sum(colors[edges[,1]] == colors[edges[,2]])
}
