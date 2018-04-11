library(boot)

# load("../raw/GCP_all_litpars/reproduceV3.Rda")
# load("../raw/GCP_PSO_iracepars/reproduceV3.Rda")
load("../raw/reproduceV3.Rda")

## helper function for bootstrapped mean
Bmean <- function(data,indices) {
  return(mean(data[indices]))
}

## Summarizing results
d <- sapply(seq_along(results.list), 
            FUN = function(i) { 
              x <- results.list[[i]]
              
              solved <- sum(x$violation == 0)
              if (solved < 100 && solved > 0) {
                # do not include evaluation cost of graphs that were not solved
                evals <- x$evaluation[x$evaluation < max(x$evaluation)]
              }
              else {
                evals <- x$evaluation
              }
              
              eval.mean <- mean(evals)
              eval.boot <- boot(data = evals, statistic=Bmean, R=1000)
              
              if (length(unique(evals)) > 1) {
                eval.ci <- boot.ci(eval.boot, type = "norm")$normal
              } else {
                eval.ci <- c(eval.mean, eval.mean, eval.mean)
              }
            c(solved = solved/100, evaluation = eval.mean, evaluation.ci.low = eval.ci[2], evaluation.ci.high = eval.ci[3])
            }
            )

d <- t(d)
d <- data.frame(d)

## Adding columns for the experiment parameters

graph.type <- rep("unset", length(results.list))
graph.type[grep("minton", names(results.list))] <- "minton"
graph.type[grep("random", names(results.list))] <- "random"

graph.size <- rep(0, length(results.list))
graph.size[grep("90", names(results.list))] <- 90
graph.size[grep("120", names(results.list))] <- 120
graph.size[grep("150", names(results.list))] <- 150
graph.size[grep("180", names(results.list))] <- 180

graph.d <- rep(0, length(results.list))
graph.d[grep("20$", names(results.list))] <- 2
graph.d[grep("25$", names(results.list))] <- 2.5
graph.d[grep("30$", names(results.list))] <- 3

method <- rep("unset", length(results.list))
method[grep("kui", names(results.list))] <- "D-ABC"
# method[grep("togashi", names(results.list))] <- "ABC.t"
method[grep("cs", names(results.list))] <- "D-CS"
method[grep("hdpso_", names(results.list))] <- "HDPSO"
method[grep("hdpso.irace", names(results.list))] <- "Tuned HDPSO"
method[grep("abc.irace", names(results.list))] <- "Tuned D-ABC"
method[grep("cs.irace", names(results.list))] <- "Tuned D-CS"
method[grep("ffa.lit", names(results.list))] <- "M-FFA"

## Adding the columns together, and sorting by experiment parameter and then method

d <- cbind(d, graph.type, graph.size, graph.d, method)
d <- d[d$method != "unset",]
d <- d[rev(with(d, order(graph.type, -graph.size, -graph.d, method))), ]


palette(c("black","red","green","blue"))
adjustcolor(7, red.f = 0)

## Some plots and analysis
# plot(d$solved)
# abline(v = 1:23*4+0.5)
# plotCI(1:96, d$evaluation, ui = d$evaluation.ci.high, li = d$evaluation.ci.low)
# abline(v = 1:23*4+0.5)


plot.success <- function(d, node) {
  d <- droplevels(d)
  
  dd <- subset(d, graph.size == node)
  for (i in levels(d$method)) {
    idx <- order(levels(d$method))[levels(d$method) == i]
  
    # ddd <- subset(dd, method == i & graph.type == "random")
    # plot(ddd$graph.d, ddd$solved,  
    #      xaxt = "n", yaxt = "n", xlab = "", ylab = "",
    #      type = "o", col = idx, pch = idx, lty = 1,
    #      ylim = c(0,1), xlim = c(2,3))
    # par(new = T)

    ddd <- subset(dd, method == i & graph.type == "minton")
    plot(ddd$graph.d, ddd$solved, 
         xaxt = "n", yaxt = "n", xlab = "", ylab = "",
         type = "o", col = idx, pch = idx, lty = 2,
         ylim = c(0,1), xlim = c(2,3))
    par(new = T)
  }
  axis(1, at = c(2,2.5,3)); axis(2, at = c(0:5) / 5)

  title(main = paste("Success Rate for N =",node,sep = " "), xlab = "Graph Density", ylab = "Solved Proportion")
  
  legend("left", bty = "n", levels(d$method), cex = 0.8,
         col = seq_along(levels(d$method)), 
         pch = seq_along(levels(d$method)))  
  # legend("bottomright", c(levels(d$method),"","Random Graph","Minton Graph"), cex = 0.8,
  #        col = c(seq_along(levels(d$method)),NA,1,1), 
  #        pch = c(seq_along(levels(d$method)),NA,NA,NA), 
  #        lty = c(seq_along(levels(d$method))^0,NA,1,2))

  par(new = F)
}

plot.evaluation <- function(d, nodes) {
  d <- droplevels(d)
  
  dd <- subset(d, graph.size == nodes)
  
  plot(0,0,type="n", main = paste("Total Evaluations for N =",nodes),
       xlim=c(1.8,3.2),ylim=c(0,max(dd$evaluation.ci.high)),
       ylab="", cex.axis = 0.8, 
       xlab="Density",
       xaxt = "n", 
       las=1)

  for (i in levels(dd$method)) {
    idx <- order(levels(dd$method))[levels(dd$method) == i]
    ddd <- subset(dd, method == i & graph.type == "minton")
    for (j in seq_len(nrow(ddd))){
      points(c(ddd$graph.d[j] - 0.04 + idx / 50, ddd$graph.d[j] - 0.04 + idx / 50), 
             c(ddd$evaluation.ci.low[j], ddd$evaluation.ci.high[j]),
             type = "l", col = idx, lty = 1)
      points(ddd$graph.d[j] - 0.04 + idx / 50, 
             ddd$evaluation[j], col = idx, pch = idx, lty = 1)
      points(c(ddd$graph.d[j] - 0.04 + idx / 50 - 0.005, ddd$graph.d[j] - 0.04 + idx / 50 + 0.005),
             c(ddd$evaluation.ci.low[j], ddd$evaluation.ci.low[j]),
             type = "l", col = idx)
      points(c(ddd$graph.d[j] - 0.04 + idx / 50 - 0.005, ddd$graph.d[j] - 0.04 + idx / 50 + 0.005),
             c(ddd$evaluation.ci.high[j], ddd$evaluation.ci.high[j]),
             type = "l", col = idx)

      }

    # ddd <- subset(dd, method == i & graph.type == "random")
    # for (j in seq_len(nrow(ddd))){
    #   points(c(ddd$graph.d[j] - 0.033 + idx / 50, ddd$graph.d[j] - 0.033 + idx / 50), 
    #          c(ddd$evaluation.ci.low[j], ddd$evaluation.ci.high[j]),
    #          type = "l", col = idx, lty = 2)
    #   points(c(ddd$graph.d[j] - 0.033 + idx / 50 - 0.005, ddd$graph.d[j] - 0.033 + idx / 50 + 0.005),
    #          c(ddd$evaluation.ci.low[j], ddd$evaluation.ci.low[j]),
    #          type = "l", col = idx)
    #   points(c(ddd$graph.d[j] - 0.033 + idx / 50 - 0.005, ddd$graph.d[j] - 0.033 + idx / 50 + 0.005),
    #          c(ddd$evaluation.ci.high[j], ddd$evaluation.ci.high[j]),
    #          type = "l", col = idx)
    # }
  }

  axis(1, at = c(2,2.5,3));

  # legend("topleft", levels(d$method), cex = 0.8,
  #        col = seq_along(levels(d$method)), 
  #        pch = seq_along(levels(d$method)))
  legend("left", bty = "n", levels(d$method), cex = 0.8,
         col = seq_along(levels(d$method)), 
         pch = seq_along(levels(d$method)))
  # legend("topleft", c(levels(d$method),"","Random Graph","Minton Graph"), cex = 0.8,
  #        col = c(seq_along(levels(d$method)),NA,1,1), 
  #        pch = c(seq_along(levels(d$method)),NA,NA,NA), 
  #        lty = c(seq_along(levels(d$method))^0,NA,1,2))
}


