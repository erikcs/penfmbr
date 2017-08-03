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
#' `vcov` can return either Shanken or GMM standard errors. For a `penfmb`
#' instance, the bootstrapped shrinkage rate is returned.
#'
#' @docType package
#' @name penfmbr
#'@importFrom Rcpp evalCpp
#'@useDynLib penfmbr
NULL


#' FMB allowing for an unbalanced panel
#' @export
#' @param factors The factors [factor1, factor2, ..., factork]
#' @param returns The returns [ret1, ret2, ..., retn]
#' @param intercept, default TRUE: fit a constant in the second stage
#' @return An `fmb` instance
fmb = function(factors, returns, intercept=TRUE) {
  factors = as.matrix(factors)
  returns = as.matrix(returns)
  stopifnot(nrow(factors) == nrow(returns))
  avg_ret = colMeans(returns, na.rm = TRUE)

  # first_stage  = .lm.fit(cbind(1, factors), returns)
  # For allowing missing values in returns
  first_stage = list()
  k = ncol(factors)
  n = ncol(returns)
  len = nrow(returns)
  coefficients = matrix(NA, k+1, n)
  residuals = matrix(NA, len, n)
  for (i in 1:n) {
    fit = lm(i ~ ., data.frame(factors, i=returns[, i]))
    coefficients[, i] = fit$coefficients
    residuals[, i] = t(returns[,i] - as.matrix(cbind(1, factors)) %*% coefficients[, i])
  }
  first_stage$coefficients = coefficients
  first_stage$residuals = residuals
  betas = as.matrix(t(first_stage$coefficients)[, -1])
  colnames(betas)= colnames(factors)

  if  (intercept) {
    # Just a dummy, allow for unbalanced panel
    fit = lm(avg_ret ~ ., data.frame(betas, avg_ret))
    lambdas = matrix(, len, k + 1)
    colnames(lambdas) = names(fit$coefficients)
    for (i in 1:len) {
      mask = complete.cases(returns[i, ])
      lambdas[i, ] = .lm.fit(cbind(1, betas[mask, ]), returns[i, mask])$coefficients
    }
    fit$coefficients = colMeans(lambdas)
    fit$fitted.values = cbind(1, betas) %*% fit$coefficients
    fit$residuals = avg_ret - fit$fitted.values
  }
  else {
    fit = lm(avg_ret ~ . -1, data.frame(betas, avg_ret))
    lambdas = matrix(, len, k)
    colnames(lambdas) = names(fit$coefficients)
    for (i in 1:len) {
      mask = complete.cases(returns[i, ])
      lambdas[i, ] = .lm.fit(betas[mask, ], returns[i, mask])$coefficients
    }
    fit$coefficients = colMeans(lambdas)
    fit$fitted.values = betas %*% fit$coefficients
    fit$residuals = avg_ret - fit$fitted.values
  }

  fit$first_stage = first_stage
  fit$avg_ret = avg_ret
  fit$betas = betas
  fit$factors = factors
  fit$returns = returns
  fit$predicted = predict(fit)
  fit$portfolio_names = colnames(returns)

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
#' @export
vcov.fmb = function(object, type = 'Shanken', ...) {
  if (type == 'Shanken')
    return (vcov_shanken(object, ...))
  if (type == 'GMM')
    return (vcov_gmm(object, ...))
  else
    stop('Unknown covariance type')
}


# Heavily overloaded to easily produce table output
#' @export
vcov.penfmb = function(object, ...) {
  shrinkage_rate = rowSums(object$bootcoefs == 0) / ncol(object$bootcoefs)
  shrinkage_rate = diag(shrinkage_rate^2)
  colnames(shrinkage_rate) = rownames(shrinkage_rate) = names(object$coefficients)

  shrinkage_rate
}


vcov_shanken = function(object, ...) {
# Cochrane (2005)
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
    cov(object$first_stage$residuals, use='complete.obs')) *
      scaling +
      cbind(0, rbind(0, var(factors)))) /
    nrow(factors)
  rownames(vcovShanken)[1] = '(Intercept)'
  colnames(vcovShanken)[1] = '(Intercept)'
  if (!intercept)
    vcovShanken = as.matrix(vcovShanken[-1, -1])

  vcovShanken
}


vcov_gmm = function(object, ...) {
# Cochrane (2005)
  covna = function(x) {
  complete = complete.cases(x)
  as.matrix(t(x[complete,])) %*% as.matrix(x[complete,]) /
    sum(complete)
  }

  ut = object$first_stage$residuals
  nT = sum(complete.cases(ut))
  factors = object$factors
  returns = object$returns
  k = ncol(factors)
  n = ncol(returns)
  # Eff = cov(factors, use='complete.obs')
  Eff = covna(factors)
  lambdas = object$coefficients
  s = k
  if (attr(terms(object), 'intercept') == 0) {
    lambdas = c(0, object$coefficients)
    s = k - 1
  }

  for (i in 1:k)
    ut = cbind(ut, object$first_stage$residuals * factors[, i])
  ut = cbind(ut, as.data.frame(returns) - t(object$first_stage$coefficients) %*% lambdas)
  ut_demeaned = ut - colMeans(ut, na.rm = TRUE)
  # S = cov(ut_demeaned, use='complete.obs')
  S = covna(ut_demeaned)
  upper = cbind(diag(n*(k+1)), matrix(0, n*(k+1), n))
  lower = cbind(matrix(0, k+1, n*(k+1)), rbind(1, t(object$betas)))
  a = rbind(upper, lower)
  d = rbind(cbind(1, t(colMeans(factors))), cbind(colMeans(factors), Eff))
  d = kronecker(d, diag(n))
  upper = cbind(d, matrix(0, n*(k+1), k+1))
  lower = cbind(matrix(0, n, n), t(kronecker(lambdas[-1], diag(n))),  cbind(1, object$betas))
  d = -rbind(upper, lower)
  vargmm = MASS::ginv(a%*%d) %*% a %*% S %*% t(a)%*%t(MASS::ginv(a%*%d)) / nT
  vcovgmm = vargmm[(ncol(vargmm)-s):ncol(vargmm), (ncol(vargmm)-s):ncol(vargmm)]
  colnames(vcovgmm) = names(object$coefficients)
  rownames(vcovgmm) = names(object$coefficients)

  vcovgmm
}


#' chisq_fmb
#' @export
#' @param object The fitted fmb model
#' @return A list containing the chi square test statistic and the p value
#' of all pricing errors jointly being zero
chisq_fmb = function(object) {
  factors = object$factors
  intercept = TRUE
  lambdas = object$coefficients[-1]
  if (attr(terms(object), 'intercept') == 0) {
    intercept = FALSE
    lambdas = object$coefficients
  }
  .sandwich = function(bread, meat)
    t(bread) %*% meat %*% bread
  scaling = 1 + as.numeric(.sandwich(lambdas, MASS::ginv(cov(factors, use='complete.obs'))))
  N = length(object$avg_ret)
  K = ncol(factors)
  alphas = matrix(residuals(object), nrow = N, ncol = 1)
  sigma = cov(object$first_stage$residuals, use='complete.obs')
  sigma_f =  cbind(0, rbind(0, var(factors))) / nrow(factors)
  I_N = diag(N)
  betas = object$betas
  betas_a = cbind(1, betas)
  left = I_N - betas_a %*% MASS::ginv( t(betas_a) %*% betas_a ) %*% t(betas_a)
  right = t(left)
  Cov_a = left %*% sigma %*% right * scaling / nrow(object$returns)
  chi = t(alphas) %*% MASS::ginv(Cov_a) %*% alphas
  pval = pchisq(chi, N - K, lower.tail = FALSE)

  list(statistic = chi, p.value = pval, mape = mean(abs(alphas)))
}
