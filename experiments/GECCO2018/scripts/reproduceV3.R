# Script for reproducing the results from Kanoh Lab in the Graph 3-GCP experiments

library(devtools)
library(EvoGCP)
library(parallel)

# Defining the graphs to be used

graphpar.list = list()
# Random Graphs
graphpar.list[["random90.20"]] = list(name = "random90.20", size = 100, nodes = 90, density = 2, method = "random")
graphpar.list[["random90.25"]] = list(name = "random90.25", size = 100, nodes = 90, density = 2.5, method = "random")
graphpar.list[["random90.30"]] = list(name = "random90.30", size = 100, nodes = 90, density = 3, method = "random")

graphpar.list[["random120.20"]] = list(name = "random120.20", size = 100, nodes = 120, density = 2, method = "random")
graphpar.list[["random120.25"]] = list(name = "random120.25", size = 100, nodes = 120, density = 2.5, method = "random")
graphpar.list[["random120.30"]] = list(name = "random120.30", size = 100, nodes = 120, density = 3, method = "random")

graphpar.list[["random150.20"]] = list(name = "random150.20", size = 100, nodes = 150, density = 2, method = "random")
graphpar.list[["random150.25"]] = list(name = "random150.25", size = 100, nodes = 150, density = 2.5, method = "random")
graphpar.list[["random150.30"]] = list(name = "random150.30", size = 100, nodes = 150, density = 3, method = "random")

graphpar.list[["random180.20"]] = list(name = "random180.20", size = 100, nodes = 180, density = 2, method = "random")
graphpar.list[["random180.25"]] = list(name = "random180.25", size = 100, nodes = 180, density = 2.5, method = "random")
graphpar.list[["random180.30"]] = list(name = "random180.30", size = 100, nodes = 180, density = 3, method = "random")

# Minton Graphs
graphpar.list[["minton90.20"]] = list(name = "minton90.20", size = 100, nodes = 90, density = 2, method = "minton")
graphpar.list[["minton90.25"]] = list(name = "minton90.25", size = 100, nodes = 90, density = 2.5, method = "minton")
graphpar.list[["minton90.30"]] = list(name = "minton90.30", size = 100, nodes = 90, density = 3, method = "minton")

graphpar.list[["minton120.20"]] = list(name = "minton120.20", size = 100, nodes = 120, density = 2, method = "minton")
graphpar.list[["minton120.25"]] = list(name = "minton120.25", size = 100, nodes = 120, density = 2.5, method = "minton")
graphpar.list[["minton120.30"]] = list(name = "minton120.30", size = 100, nodes = 120, density = 3, method = "minton")

graphpar.list[["minton150.20"]] = list(name = "minton150.20", size = 100, nodes = 150, density = 2, method = "minton")
graphpar.list[["minton150.25"]] = list(name = "minton150.25", size = 100, nodes = 150, density = 2.5, method = "minton")
graphpar.list[["minton150.30"]] = list(name = "minton150.30", size = 100, nodes = 150, density = 3, method = "minton")

graphpar.list[["minton180.20"]] = list(name = "minton180.20", size = 100, nodes = 180, density = 2, method = "minton")
graphpar.list[["minton180.25"]] = list(name = "minton180.25", size = 100, nodes = 180, density = 2.5, method = "minton")
graphpar.list[["minton180.30"]] = list(name = "minton180.30", size = 100, nodes = 180, density = 3, method = "minton")

# Defining Methods
solverpar.list = list()
solverpar.list[["cs.toda"]] = list(name = "cs.toda", method = "cuckoo", pop = 10, pc = 0.0001, compare = F, policy = "levy", alpha = 1, beta = 1.5, E = NULL)
solverpar.list[["abc.togashi"]] = list(name = "abc.togashi", method = "abc", pop = 500, onlooker = 500, c = 3, scout = 1, limit = 300)
solverpar.list[["hdpso"]] = list(name = "hdpso", method = "hdpso", pop = 10, w = 0.05, c1 = 7, c2 = 0.03)
solverpar.list[["abc.kui"]] = list(name = "abc.kui", method = "abc_kui", pop = 200, onlooker = 200, c = 2, scout = 1, limit = 90)

tester.batch(solverpar.list, graphpar.list, results.file = "reproduceV3.Rda", graphs.file = "graphs.90.to.180.Rda", nfe = 20000000, random.seed = 42, ncores = 12)
