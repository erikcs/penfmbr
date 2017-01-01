library(penfmbr)

factors = read.csv('factors.csv')[, -1]
returns = read.csv('returns.csv')[, -1]

test_that("fmb coef and vcov", {
  out = fmb(factors[, 1, drop = FALSE], returns)
  expected = c(1.4744, -0.6972)
  expect_equal(out$coefficients, expected, check.names = FALSE, tolerance = 1e-3)
  expect_equal(sqrt(diag(vcov(out))), c(0.4422, 0.4832),
               check.names = FALSE, tolerance = 1e-3)

  out = fmb(factors[, c(1, 2, 3)], returns)
  expected = c(1.351, -0.7945, 0.1415, 0.4257)
  expect_equal(out$coefficients, expected, check.names = FALSE, tolerance = 1e-3)
  expect_equal(sqrt(diag(vcov(out))), c(0.2898, 0.3545, 0.1428, 0.1354),
               check.names = FALSE, tolerance = 1e-3)
})

test_that("penfmb", {
  out = penfmb(factors[, 1, drop = FALSE], returns, nboot = 2)
  expected = c(1.4744, -0.6972)
  expect_equal(out$coefficients, expected, check.names = FALSE, tolerance = 1e-3)

  out = penfmb(factors, returns, nboot = 2)
  expected = c(1.2988, -0.7470, 0.1425, 0.4294, 0.0)
  expect_equal(out$coefficients, expected, check.names = FALSE, tolerance = 1e-3)
})

test_that("coordinate descent equals lm", {
  y = as.matrix(mtcars)[, 'mpg']
  X = as.matrix(mtcars)[, c('cyl', 'disp', 'hp')]
  out = cd(X, y, 0., rep(1, ncol(X)), 1e-6, 1e5)
  expected = coef(lm(y ~ X))

  expect_equal(unlist(out), expected, check.names = FALSE, tolerance = 1e-6)
})

test_that("coordinate descent equals glmnet", {
  y = as.matrix(mtcars)[, 'mpg']
  X = as.matrix(mtcars)[, c('cyl', 'disp', 'hp')]
  out = cd(X, y, 0.5, rep(1, ncol(X)), 1e-6, 1e5)
  expected = c(32.8, -1.16, -0.017, -0.01)

  expect_equal(unlist(out), expected, check.names = FALSE, tolerance = 1e-2)
})
