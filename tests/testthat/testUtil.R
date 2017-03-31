library(EvoGCP)
context("Utility Functions")

set.seed(42)

test_that("Generate Random Coloring", {
  c = random_coloring(10)
  expect_equal(length(c), 10)
  expect_lte(max(c),3)
  expect_gte(min(c),1)
})

test_that("Graph Connectivity is correct", {
  V <- 5
  E1 <- matrix(c(1, 2, 1, 3, 1, 4, 1, 5), 4, byrow = T)
  E2 <- matrix(c(1, 2, 1, 3), 2, byrow = T)
  expect_true(connectivity(V, E1))
  expect_false(connectivity(V, E2))
})

test_that("Levy Distribution",
{
  rl1 <- rlevy (1, 1)
  rl1000 <- rlevy (1000, 1)
  expect_type(rl1, "double")
  expect_length(rl1, 1)
  expect_type(rl1000, "double")
  expect_length(rl1000, 1000)
})

test_that("Similarity", {
  x <- c(1,2,2,1,3)
  y <- c(1,1,1,1,1)
  xmat <- matrix(c(x, y), 4, 5, byrow = T) # x, y, x, y
  ymat <- matrix(c(x, x), 4, 5, byrow = T) # x, x, x, x
  expect_equal(similarity(x, x), 1)
  expect_equal(similarity(x, y), 0.4)
  expect_equal(similarity(xmat, ymat), c(1, 0.4, 1, 0.4))
  expect_equal(similarity(xmat, xmat), c(1, 1, 1, 1))
})

test_that("Kui's Similarity", {
  x <- c(1,2,2,1,3)
  y <- c(1,1,1,1,1)
  xmat <- matrix(c(x, y), 4, 5, byrow = T) # x, y, x, y
  ymat <- matrix(c(x, x), 4, 5, byrow = T) # x, x, x, x
  # Kui's similarity adds a random factor)
  set.seed(42)
  expect_equal(similarity.kui(x, x), 1)
  expect_lte(similarity.kui(x, y) - 0.438, 0.001)
  expect_lte(sum(similarity.kui(xmat, ymat) - c(1, 0.502, 1, 0.689)), 0.001)
  expect_equal(similarity.kui(xmat, xmat), c(1, 1, 1, 1))
})


