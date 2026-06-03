test_that("get_mcc returns 1 for perfect recovery", {
  expect_equal(get_mcc(1:5, 1:5, p = 20), 1)
})

test_that("get_mcc returns 0 for empty selection", {
  expect_equal(get_mcc(1:5, integer(0), p = 20), 0)
})

test_that("get_laplacian returns correct Laplacian", {
  A <- matrix(c(0,1,0, 1,0,1, 0,1,0), 3, 3)
  L <- get_laplacian(A)
  expect_equal(rowSums(L), c(0, 0, 0))
})

test_that("make_graph returns symmetric matrix", {
  A <- make_graph(p = 20, type = "chain")
  expect_equal(A, t(A))
})

test_that("generate_data returns correct dimensions", {
  set.seed(1)
  A    <- make_graph(p = 10, type = "chain")
  bobj <- generate_beta(A, s2 = 3, signal = 2)
  dat  <- generate_data(n = 50, p = 10, q = 30,
                        beta_true = bobj$beta_true)
  expect_equal(nrow(dat$X), 50)
  expect_equal(ncol(dat$X), 10)
  expect_equal(length(dat$Y), 50)
})
