library(EvoGCP)
context("Generating 3-GCP Problems")

set.seed(42)

test_that("Generate Random gives valid graph", {
  g = graph_random(10, 2)
  expect_equal(dim(g$E), c(round(10 * 2), 2))
  expect_equal(sum(duplicated(g$E)), 0)
  expect_lte(max(g$E),10)
  expect_gte(min(g$E),0)
})

test_that("Generate Minton gives valid graph", {
  g = graph_minton(10, 2)
  expect_equal(dim(g$E), c(round(10 * 2), 2))
  expect_equal(sum(duplicated(g$E)), 0)
  expect_lte(max(g$E),10)
  expect_gte(min(g$E),0)
})

test_that("Generate Culberson Flat gives valid graph", {
  g = graph_culberson_flat(10,0.014,0)
  expect_equal(ncol(g$E), 2)
  expect_equal(sum(duplicated(g$E)), 0)
  expect_lte(max(g$E),10)
  expect_gte(min(g$E),0)

  #test coloring
  g$adj <- adjacency_list(g)
  degree <- sapply(1:g$V, FUN = function(x) { sum(g$E[1]==x)+sum(g$E[2]==x) })
  s <- solver_dsatur(G, G$V+1, args=list(permutation=degree))
  expect_equal(s$violation, 0)
  expect_equal(length(s$best[s$best==0]), 0)
})

test_that("Generate Culberson Equipartite gives valid graph", {
  g = graph_culberson_equipartite(10, .5)
  expect_equal(ncol(g$E), 2)
  expect_equal(sum(duplicated(g$E)), 0)
  expect_lte(max(g$E),10)
  expect_gte(min(g$E),0)

  #test coloring
  g$adj <- adjacency_list(g)
  degree <- sapply(1:g$V, FUN = function(x) { sum(g$E[1]==x)+sum(g$E[2]==x) })
  s <- solver_dsatur(G, G$V+1, args=list(permutation=degree))
  expect_equal(s$violation, 0)
  expect_equal(length(s$best[s$best==0]), 0)
})

test_that("Generate Culberson with Variability valid graph", {
  g = graph_culberson_variability(10, .3, 2)
  expect_equal(ncol(g$E), 2)
  expect_equal(sum(duplicated(g$E)), 0)
  expect_lte(max(g$E),10)
  expect_gte(min(g$E),0)

  #test coloring
  g$adj <- adjacency_list(g)
  degree <- sapply(1:g$V, FUN = function(x) { sum(g$E[1]==x)+sum(g$E[2]==x) })
  s <- solver_dsatur(G, G$V+1, args=list(permutation=degree))
  expect_equal(s$violation, 0)
  expect_equal(length(s$best[s$best==0]), 0)
})

test_that("Generate Sets of Graphs", {
  ps = problemset(10, 100, 2.5, "random")
  expect_equal(length(ps), 10)
  ps = problemset(10, 100, 2.5, "minton")
  expect_equal(length(ps), 10)
  ps = problemset(10, 100, 1, "culberson_equipartite")
  expect_equal(length(ps), 10)
  ps = problemset(10, 100, 1, x=0, "culberson_flat")
  expect_equal(length(ps), 10)
  ps = problemset(10, 100, 1, x=0, "culberson_variability")
  expect_equal(length(ps), 10)
})
