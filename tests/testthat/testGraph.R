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

test_that("Generate Sets of Graphs", {
  ps = problemset(10, 100, 2.5, "random")
  expect_equal(length(ps), 10)
  ps = problemset(10, 100, 2.5, "minton")
  expect_equal(length(ps), 10)
})
