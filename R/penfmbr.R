#' penfmbr: Plain FamaMacBeth and Penalized FamaMacBeth estimators.
#'
#' For the penalized FamaMacBeth, the level tuning parameter is set to the
#' average of the residual standard deviation from the first stage. The
#' dimensionalty paramter is 4 by default, and individual penalty weights are
#' calculated as the Fisher transform of partial correlations.
#'
#' `returns` and `factors` are two containers with equal number of observations
#' (aligned) that can be coerced to a matrix.
#' Both estimators return a hacked together sublass of `lm` that plays
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
fmb = function(factors, returns) {
  factors = as.matrix(factors)
  returns = as.matrix(returns)
  avg_ret = colMeans(returns)

  first_stage  = .lm.fit(cbind(1, factors), returns)
  betas = as.matrix(t(first_stage$coefficients)[, -1])
  colnames(betas)= colnames(factors)

  fit = lm(avg_ret ~ ., data.frame(betas, avg_ret))
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
#' @param alpha Tthe penalty parameter, default 4
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
    fit = cd(betas, avg_ret, alpha = nu, alpha_weights = penalty, tol=1e-6, maxiter = 1e6)

    c(fit$w0, fit$w)
  }

  factors = as.matrix(factors)
  returns = as.matrix(returns)

  # Get the shrinkage rate
  bootidx = tseries::tsbootstrap(1:nrow(returns), nb = nboot, b = nblocks, type = 'stationary')
  bootcoefs = apply(bootidx, 2, function(x)
    .penfmb.fit(factors[x, ,drop = FALSE], returns[x, ], alpha))

  # Dummy obj
  dfit = fmb(factors, returns)
  coefs = .penfmb.fit(factors, returns, alpha)

  dfit$coefficients = coefs
  dfit$fitted.values = cbind(1, dfit$betas) %*% dfit$coefficients
  dfit$residuals = dfit$avg_ret - dfit$fitted.values
  dfit$bootcoefs = bootcoefs

  names(dfit$coefficients) = c('(Intercept)', colnames(factors))
  class(dfit)[1] = 'penfmb'

  dfit
}


# Shanken standard errors
# TODO: vcov from GMM estimation
#' @export
vcov.fmb = function(object, ...) {
  .sandwich = function(bread, meat)
    t(bread) %*% meat %*% bread

  betas = object$betas
  lambdas = object$coefficients[-1]
  factors = object$factors

  scaling = 1 + as.numeric(.sandwich(lambdas, MASS::ginv(var(factors))))

  vcovShanken = (.sandwich(
    MASS::ginv(cbind(1, betas) %*% t(cbind(1, betas))) %*% cbind(1, betas),
    var(object$first_stage$residuals)) *
      scaling +
      cbind(0, rbind(0, var(factors)))) /
    nrow(factors)

  rownames(vcovShanken)[1] = '(Intercept)'
  colnames(vcovShanken)[1] = '(Intercept)'

  vcovShanken
}


# Heavily overloaded
#' @export
vcov.penfmb = function(object, ...) {

  shrinkage_rate = rowSums(object$bootcoefs == 0) / ncol(object$bootcoefs)
  vcov = diag(shrinkage_rate^2)
  colnames(vcov) = rownames(vcov) = names(object$coefficients)

  vcov
}
