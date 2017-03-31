library(EvoGCP)
context("3-GCP Problem works correctly")

set.seed(42)

edges = list(N = 4, E = t(matrix(c(1, 2, 1, 3, 2, 3), 2)))
colors1 = c(1,3,2,1)
colors2 = c(1,1,1,1)

test_that("Fitness Function returns correct values for simple cases", {
  expect_equal(evaluate(colors1, edges), 0)
  expect_equal(evaluate(colors2, edges), 3)
})

test_that("Evaluation evaluates correctly a generated case", {
  expect_gte(evaluate(colors1, graph_random(4, 1)), 0)
})
