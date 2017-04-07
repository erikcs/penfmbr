#' penfmbr: Plain FamaMacBeth and Penalized FamaMacBeth estimators.
#'
#' For the penalized FamaMacBeth, the level tuning parameter is set to the
#' average of the residual standard deviation from the first stage. The
#' dimensionalty paramter is 4 by default, and individual penalty weights are
#' calculated as the Fisher transform of partial correlations.
#'
#' `returns` and `factors` are two containers with equal number of observations
#' (aligned) that can be coerced to a matrix.
#' Both estimators return a crudely hacked together sublass of `lm` that plays
#' nicely with `stargazer`.
#'
#' @docType package
#' @name penfmbr
#'@importFrom Rcpp evalCpp
#'@useDynLib penfmbr
NULL

#' fmb
#' @export
#' @param factors The factors
#' @param returns The returns
#' @return An `fmb` instance
fmb = function(factors, returns, intercept=TRUE) {
  factors = as.matrix(factors)
  returns = as.matrix(returns)
  avg_ret = colMeans(returns)

  first_stage  = .lm.fit(cbind(1, factors), returns)
  betas = as.matrix(t(first_stage$coefficients)[, -1])
  colnames(betas)= colnames(factors)

  if (intercept)
    fit = lm(avg_ret ~ ., data.frame(betas, avg_ret))
  else
    fit = lm(avg_ret ~ . -1, data.frame(betas, avg_ret))

  fit$first_stage = first_stage
  fit$avg_ret = avg_ret
  fit$betas = betas
  fit$factors = factors
  fit$returns = returns

  class(fit) = c('fmb', 'lm')
  fit
}

#' penfmb
#' @export
#' @param factors The factors
#' @param returns The returns
#' @param alpha The penalty parameter, default 4
#' @param nboot Number of bootstrap replications for the shrinkage rate, default 100
#' @param nblocks Number of blocks in the stationary bootstrap, default 5
#' @return A `penfmb` instance
penfmb = function(factors, returns, alpha = 4, nboot=100, nblocks=5) {

  .penalty = function(factors, returns, alpha) {
    # Fisher transform of partial correlations
    pc = t(apply(returns, 2, function(x)
      -cov2cor(solve(cov(cbind(x, factors))))[-1, 1]))

    if (ncol(pc) != ncol(factors))
      pc = t(pc)

    colSums(abs(0.5 * log((1 + pc) / (1 - pc)))) ^ (-alpha) *
      ncol(returns)^alpha / nrow(returns)^(alpha/2)
  }

  .penfmb.fit = function(factors, returns, alpha) {

    avg_ret = colMeans(returns)
    first_stage  = .lm.fit(cbind(1, factors), returns)
    betas = as.matrix(t(first_stage$coefficients)[, -1])
    nu = sqrt(mean(diag(var(first_stage$residuals))))
    penalty = .penalty(factors, returns, alpha)

    # No good reason to avoid using `glmnet` except that it did not play nicely
    # with all penalties = 0, and 1-column `X` input.
    fit = cd(betas, avg_ret, alpha = nu, alpha_weights = penalty, tol=1e-5, maxiter = 1e5)

    c(fit$w0, fit$w)
    # fit = glmnet::glmnet(betas, avg_ret, alpha = 1, lambda = nu, penalty.factor = penalty)
    # c(as.vector(fit$a0), as.vector(fit$beta))
  }

  factors = as.matrix(factors)
  returns = as.matrix(returns)

  # Get the shrinkage rate
  bootidx = tseries::tsbootstrap(1:nrow(returns), nb = nboot, b = nblocks, type = 'stationary')
  bootcoefs = apply(bootidx, 2, function(x)
    .penfmb.fit(factors[x, , drop = FALSE], returns[x, ], alpha))

  # Dummy obj
  dfit = fmb(factors, returns)
  dfit$coefficients = .penfmb.fit(factors, returns, alpha)

  dfit$fitted.values = cbind(1, dfit$betas) %*% dfit$coefficients
  dfit$residuals = dfit$avg_ret - dfit$fitted.values
  dfit$bootcoefs = bootcoefs

  names(dfit$coefficients) = c('(Intercept)', colnames(factors))
  class(dfit)[1] = 'penfmb'

  dfit
}


# Shanken standard errors by default
# TODO: maybe GMM vcov
#' @export
vcov.fmb = function(object, ...) {
  .sandwich = function(bread, meat)
    t(bread) %*% meat %*% bread

  intercept = TRUE
  lambdas = object$coefficients[-1]
  if (attr(terms(object), 'intercept') == 0) {
    intercept = FALSE
    lambdas = object$coefficients
  }

  betas = object$betas
  factors = object$factors

  scaling = 1 + as.numeric(.sandwich(lambdas, MASS::ginv(cov(factors))))

  vcovShanken = (.sandwich(
    MASS::ginv(cbind(1, betas) %*% t(cbind(1, betas))) %*% cbind(1, betas),
    cov(object$first_stage$residuals)) *
      scaling +
      cbind(0, rbind(0, var(factors)))) /
    nrow(factors)

  rownames(vcovShanken)[1] = '(Intercept)'
  colnames(vcovShanken)[1] = '(Intercept)'

  if (!intercept)
    vcovShanken = as.matrix(vcovShanken[-1, -1])

  vcovShanken
}


# Heavily overloaded to easily produce table output
#' @export
vcov.penfmb = function(object, ...) {

  shrinkage_rate = rowSums(object$bootcoefs == 0) / ncol(object$bootcoefs)
  shrinkage_rate = diag(shrinkage_rate^2)
  colnames(shrinkage_rate) = rownames(shrinkage_rate) = names(object$coefficients)

  shrinkage_rate
}
