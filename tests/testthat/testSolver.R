library(EvoGCP)
context("Solver Components")

test_that("Testing batch solver outputs", {
  set.seed(42)
  P <- problemset(10, 20, 1, "random")
  result <- tester(P, nfe = 10, solverpar = list(name = "random"))
  testthat::expect_named(result, c('violation', 'solution', 'evaluation'))
  testthat::expect_length(result$violation, 10)
  testthat::expect_length(result$evaluation, 10)
  testthat::expect_length(result$solution, 10 * 20)
})

##################################################################

test_that("Testing Cuckoo Mutation", {
  set.seed(42)
  x <- sample(3, 50, replace = T)
  xt <- cuckoo.mutate(x, 10)
  testthat::expect_gte(min(xt), 1)
  testthat::expect_lte(max(xt), 3)
  testthat::expect_equal(sum(x != xt), 10)

  xt <- cuckoo.mutate(x, 0)
  testthat::expect_equal(x, xt)
})

test_that("Testing Cuckoo Step 1", {
  set.seed(42)
  G <- graph_random(10, 2)
  P <- t(sapply(1:5, FUN = function(x) { random_coloring(G$V) }))
  ret <- cuckoo.step1(P, G, alpha = 1, beta = 1)
  testthat::expect_named(ret, c("P", "V", "eval"))
  testthat::expect_equal(nrow(ret$P), 5)
  testthat::expect_equal(length(ret$V), 5)
  testthat::expect_equal(ret$eval, sum(ret$V < G$V+1))
})


test_that("Testing Cuckoo Step 2", {
  set.seed(42)
  G <- graph_random(10, 2)
  P <- t(sapply(1:5, FUN = function(x) { random_coloring(G$V) }))
  ret <- cuckoo.step2(P, G, alpha = 1, beta = 1, pc = 0.5, ppolicy = "fixed", E = 2)
  testthat::expect_named(ret, c("P", "V", "eval"))
  testthat::expect_equal(nrow(ret$P), 5)
  testthat::expect_equal(length(ret$V), 5)
  testthat::expect_equal(ret$eval, sum(ret$V < G$V+1))
})

####################################################

test_that("Testing ABC Mutation", {
  xi <- c(1,2,3,1,2,3)
  xj <- c(3,3,2,3,3,2)
  testthat::expect_equal(abc.mutate(xi, xi, 3), xi)
  xu <- abc.mutate(xi, xj, 4)
  testthat::expect_false(all(xi == xu))
  testthat::expect_equal(xi == xu, xj != xu)
  testthat::expect_gte(min(xu), 1)
  testthat::expect_lte(max(xu), 3)
})

test_that("Testing ABC step 2", {
  G <- list(V = 5, E = matrix(c(1, 2, 1, 3, 2, 4, 3, 5, 2, 5, 4, 5), 6, byrow = T))
  P <- matrix(rep(1:3,15),9)
  P[3,2] <- 2
  V <- apply(P, 1, evaluate, graph = G)
  set.seed(42)
  newP <- abc.step2(P, G, V, c = 3, pop.onlooker = 10)
  testthat::expect_equal_to_reference(newP, "testSolver.abcstep2.RDS")
})

test_that("Testing ABC step 3", {
  G <- graph_random(10,2)
  P <- t(sapply(1:5, FUN = function(x) { random_coloring(G$V) }))
  A <- c(0, 0, 4, 3, 5)

  set.seed(42)
  N <- abc.step3(P, G, A, limit = 10, pop.scout = 5)
  testthat::expect_equal(N$eval, 0)
  testthat::expect_equal(N$P, P)

  set.seed(42)
  N <- abc.step3(P, G, A, limit = 1, pop.scout = 1)
  testthat::expect_equal(N$eval, 1)
  testthat::expect_equal(N$P[-5, ], P[-5, ])
  testthat::expect_true(any(N$P[5, ] != P[5, ]))

  set.seed(42)
  N <- abc.step3(P, G, A, limit = 1, pop.scout = 2)
  testthat::expect_equal(N$eval, 2)
  testthat::expect_equal(N$P[c(1,2,4), ], P[c(1,2,4), ])
  testthat::expect_true(any(N$P[3, ] != P[3, ]))
  testthat::expect_true(any(N$P[5, ] != P[5, ]))

  set.seed(42)
  N <- abc.step3(P, G, A, limit = 1, pop.scout = 5)
  testthat::expect_equal(N$eval, 3)
  testthat::expect_equal(N$P[c(1,2), ], P[c(1,2), ])
  testthat::expect_true(any(N$P[3, ] != P[3, ]))
  testthat::expect_true(any(N$P[3, ] != P[4, ]))
  testthat::expect_true(any(N$P[5, ] != P[5, ]))
})

##################################################################
test_that("Testing DSatur, valid coloring", {
  set.seed(42)
  G <-  graph_random(10, 2)
  G$adj <- adjacency_list(G)
  degree <- sapply(1:G$V, FUN = function(x) { sum(G$E[,1]==x)+sum(G$E[,2]==x) })

  set.seed(42)
  s <- solver_dsatur(G, G$V+1, args=list(weight=degree, partial_solution=NULL, leave_uncolored=F))
  testthat::expect_equal(s$violation, 0)
  testthat::expect_equal(length(s$best[s$best==0]), 0)

  set.seed(42)
  s <- solver_dsatur(G, G$V+1, args=list(weight=degree, partial_solution=NULL, leave_uncolored=T))
  testthat::expect_equal(s$violation, 0)
  testthat::expect_equal(length(s$best[s$best==0]), 0)
  testthat::expect_equal(length(s$best[s$best==0]), s$violation)
})

test_that("Testing DSatur with a graph which generates violations", {
  set.seed(42)
  #G <-  graph_culberson_equipartite(500, .015)
  G <-  graph_random(45, 2)
  G$adj <- adjacency_list(G)
  degree <- sapply(1:G$V, FUN = function(x) { sum(G$E[,1]==x)+sum(G$E[,2]==x) })

  set.seed(42)
  s <- solver_dsatur(G, G$V+1, args=list(weight=degree, partial_solution=NULL, leave_uncolored=F))
  testthat::expect_equal(s$violation, 1)
  testthat::expect_equal(length(s$best[s$best==0]), 0) #violation due colored vertex

  set.seed(42)
  s <- solver_dsatur(G, G$V+1, args=list(weight=degree, partial_solution=NULL, leave_uncolored=T))
  testthat::expect_equal(s$violation, 1)
  testthat::expect_equal(length(s$best[s$best==0]), 1)
  testthat::expect_equal(length(s$best[s$best==0]), s$violation) #violation due uncolored vertex
})



##################################################################
test_that("Testing Firefly Algorithm Step 1", {
  #setup
  set.seed(42)
  G <-  graph_random(10, 2)
  G$adj <- adjacency_list(G)

  set.seed(42)
  W <- t(sapply(1:5, FUN = function(x) { random_weight(10, 0, 1) }))
  set.seed(42)
  newP <- ffa.decode(W, G)
  P <- newP$P
  V <- newP$V

  #test
  set.seed(42)
  newW <- ffa.step1(W, V, alpha=0.1, beta=0.1, gamma=0.8, nrow(W))
  testthat::expect_equal(nrow(W), nrow(newW))
  testthat::expect_equal(ncol(W), ncol(newW))
  testthat::expect_true(all(newW > 0 & newW < 1))
})

test_that("Testing Firefly Algorithm Step 2", {
  #setup
  set.seed(42)
  G <-  graph_random(10, 2)
  G$adj <- adjacency_list(G)

  set.seed(42)
  W <- t(sapply(1:5, FUN = function(x) { random_weight(10, 0, 1) }))
  set.seed(42)
  newP <- ffa.decode(W, G)
  P <- newP$P
  V <- newP$V

  #test
  set.seed(42)
  newW <- ffa.step2(P, W, G)
  testthat::expect_equal(nrow(W), nrow(newW))
  testthat::expect_equal(ncol(W), ncol(newW))
  testthat::expect_true(all(newW > 0 & newW < 1))
})

test_that("Testing Firefly Algorithm Operator Heuristical Swap", {
  #setup
  G <- list(V = 5, E = matrix(c(1, 2, 1, 3, 2, 4, 3, 5, 2, 5, 4, 5), 6, byrow = T))
  G$adj <- adjacency_list(G)
  W <- t(sapply(1:G$V, FUN = function(x) { sum(G$E[,1]==x)+sum(G$E[,2]==x) }))
  P <- t(c(1,2,3,3,0))

  #test
  newW <- ffa.step2(P, W, G)
  testthat::expect_equal(nrow(W), nrow(newW))
  testthat::expect_equal(ncol(W), ncol(newW))
  testthat::expect_true(all(newW == c(3,3,2,2,2)) | all(newW == c(2,3,2,2,3))) #it will swap 5 with 1 or 2

  #test
  set.seed(42)
  W <- t(sapply(1:5, FUN = function(x) { random_weight(5, 0, 1) }))
  set.seed(42)
  P <- ffa.decode(W, G)
  P <- P$P
  newW <- ffa.step2(P, W, G)
  testthat::expect_equal(nrow(W), nrow(newW))
  testthat::expect_equal(ncol(W), ncol(newW))
})

test_that("Testing Firefly Algorithm", {
  #setup
  set.seed(42)
  G <-  graph_random(100, 2)

  #test
  set.seed(42)
  result <- solver_ffa(G, 1000, args=list(pop=5, lb=0, ub=1, alpha=0.1, beta=0.1, gamma=0.8, swap=F))
  testthat::expect_named(result, c("violation", "best", "evals"))
  testthat::expect_length(result$violation, 1)
  testthat::expect_length(result$best, 100)
  testthat::expect_length(result$evals, 1)
  testthat::expect_gte(min(result$best), 0)
  testthat::expect_lte(max(result$best), 3)

  #swap
  set.seed(42)
  result <- solver_ffa(G, 1000, args=list(pop=5, lb=0, ub=1, alpha=0.1, beta=0.1, gamma=0.8, swap=T))
  testthat::expect_named(result, c("violation", "best", "evals"))
  testthat::expect_length(result$violation, 1)
  testthat::expect_length(result$best, 100)
  testthat::expect_length(result$evals, 1)
  testthat::expect_gte(min(result$best), 0)
  testthat::expect_lte(max(result$best), 3)

  set.seed(42)
  G <-  graph_culberson_equipartite(500, .015)

  set.seed(42)
  result <- solver_ffa(G, 30000, args=list(pop=10, lb=0, ub=1, alpha=0.1, beta=0.1, gamma=0.8, swap=F))
  testthat::expect_named(result, c("violation", "best", "evals"))
  testthat::expect_length(result$violation, 1)
  testthat::expect_length(result$best, 500)
  testthat::expect_length(result$evals, 1)
  testthat::expect_gte(min(result$best), 0)
  testthat::expect_lte(max(result$best), 3)

  set.seed(42)
  result <- solver_ffa(G, 30000, args=list(pop=10, lb=0, ub=1, alpha=0.1, beta=0.1, gamma=0.8, swap=T))
  testthat::expect_named(result, c("violation", "best", "evals"))
  testthat::expect_length(result$violation, 1)
  testthat::expect_length(result$best, 500)
  testthat::expect_length(result$evals, 1)
  testthat::expect_gte(min(result$best), 0)
  testthat::expect_lte(max(result$best), 3)

})

