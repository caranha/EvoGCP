library(EvoGCP)
context("Solver Results")

test_that("Random Solver has not changed results", {
  set.seed(42)
  G <- graph_minton(20,2)

  ret <- solver_random(G, 50, NULL)
  testthat::expect_equal_to_reference(ret, "testSolverResults.random.RDS")
})

test_that("Cuckoo Solver has not changed results", {
  set.seed(42)
  G <- graph_minton(20,2)

  r1 <- solver_cuckoo(G, 50, list(pop = 10, pc = 0.2, compare = F, policy = "uniform", alpha = 1, beta = 1, E = 5))
  r2 <- solver_cuckoo(G, 50, list(pop = 10, pc = 0.2, compare = F, policy = "fixed", alpha = 1, beta = 1, E = 5))
  r3 <- solver_cuckoo(G, 50, list(pop = 10, pc = 0.2, compare = T, policy = "fixed", alpha = 1, beta = 1, E = 5))
  r4 <- solver_cuckoo(G, 50, list(pop = 10, pc = 0.2, compare = F, policy = "levy", alpha = 1, beta = 1, E = 5))
  testthat::expect_equal_to_reference(list(r1,r2,r3,r4), "testSolverResults.cuckoo.RDS")
})

test_that("ABC Solver has not changed results", {
  set.seed(42)
  G <- graph_minton(20,2)

  ret <- solver_abc(G, 50, list(pop = 10, c = 5, onlooker = 10, scout = 2, limit = 5))
  testthat::expect_equal_to_reference(ret, "testSolverResults.abc.RDS")
})

test_that("ABC Solver (Kui) has not changed results", {
  set.seed(42)
  G <- graph_minton(20,2)

  ret <- solver_abc_kui(G, 50, list(pop = 10, c = 5, onlooker = 10, scout = 2, limit = 5))
  testthat::expect_equal_to_reference(ret, "testSolverResults.abc_kui.RDS")
})

test_that("HDPSO Solver has not changed results", {
  set.seed(42)
  G <- graph_minton(20,2)

  ret <- solver_hdpso(G, 50, list(pop = 10, w = 0.1, c1 = 2, c2 = 0.5))
  testthat::expect_equal_to_reference(ret, "testSolverResults.hdpso.RDS")
})

test_that("DSatur Solver [No Violation] has not changed results", {
  set.seed(42)
  G <- graph_minton(20,2)
  G$adj <- adjacency_list(G)
  degree <- sapply(1:G$V, FUN = function(x) { sum(G$E[,1]==x)+sum(G$E[,2]==x) })

  set.seed(42)
  ret <- solver_dsatur(G, G$V+1, args=list(weight=degree, return_satur=F))
  testthat::expect_equal_to_reference(ret, "testSolverResults.DSatur1.RDS")
})

test_that("DSatur Solver [Violation] has not changed results", {
  set.seed(42)
  G <-  graph_random(45, 2)
  G$adj <- adjacency_list(G)
  degree <- sapply(1:G$V, FUN = function(x) { sum(G$E[,1]==x)+sum(G$E[,2]==x) })

  set.seed(42)
  ret <- solver_dsatur(G, G$V+1, args=list(weight=degree, return_satur=F))
  testthat::expect_equal_to_reference(ret, "testSolverResults.DSatur2.RDS")
})

test_that("Firefly Algorithm Solver has not changed results", {
  set.seed(42)
  G <- graph_minton(20,2)
  #G <- graph_culberson_equipartite(500, .15)

  set.seed(42)
  ret <- solver_ffa(G, 50, list(pop=500, lb=0, ub=1, alpha=0.1, beta=0.1, gamma=0.8))
  testthat::expect_equal_to_reference(ret, "testSolverResults.ffa.RDS")
})
