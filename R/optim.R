###########################################################################
#
# Optimisation methods for non smooth functions
#
# Tathagata Basu & Matthias C. M. Troffaes
# 25 Sept 2019
#
###########################################################################

#' Optimisation tool for non smooth functions
#'
#' Function for minimising piecewise differentiable objective functions.
#'
#' @param theta Initial value.
#' @param fn Function to optimise.
#' @param df Subgradient of the objective function.
#' @param tol Relative convergence tolerance
#' @param t Stepsizes for updating the set of parameters.
#' @return The function returns a list with following components
#' \item{par}{The best set of parameters}
#' \item{value}{The value of the function corresponding to best set of parameters.}
#' \item{iteration}{Number of iterations taken}
#' @export

sgOptim = function(theta, fn, df, tol = 1e-6, t) {
  value = fn(theta)
  par = theta.last = theta
  maxit = length(t)
  for (i in 1:maxit) {
    theta = theta.last - t[i] * df(theta.last)
    if(is.na(sum((df(theta) ^ 2))))
      stop("the subgradient vector is returning undefined value, change stepsize")
    fx = fn(theta)
    if(fx < value) {
      par = theta
      value = fx
    }
    if(sqrt(sum((df(theta) ^ 2))) < tol)
      break()
    else
      theta.last = theta
  }

  output = list("par" = as.vector(par), "value" = value, "iteration" = i)
  output
}

#' Proximal gradient optimisation for non smooth functions
#'
#' Function for minimising piecewise differentiable objective functions.
#'
#' @param theta Initial value.
#' @param fn Function to optimise.
#' @param df Subgradient of the objective function.
#' @param prox Proximal operator for proximal gradient method.
#' @param acl Flag for accelerated proximal gradient method. Default is false.
#' @param tol Relative convergence tolerance
#' @param t Stepsizes for updating the set of parameters.
#' @return The function returns a list with following components
#' \item{par}{The best set of parameters}
#' \item{value}{The value of the function corresponding to best set of parameters.}
#' \item{iteration}{Number of iterations taken}
#' @export

pgOptim = function(theta, fn, df, prox, acl = F, tol = 1e-6, t) {
  par = theta.prev = theta
  maxit = length(t)
  if(acl){
    for (i in 1:maxit) {
      v = theta + (i - 2) * (theta - theta.prev) / (i + 1)
      par = prox(t[i], v - t[i] * df(v))
      theta.prev = theta
      theta = par
      if(is.na(sqrt(sum((df(theta) ^ 2)))))
        stop("the subgradient vector is returning undefined value, change stepsize")
      if(sqrt(sum((df(theta) ^ 2))) < tol)
        break()
    }
  }
  else{
    for (i in 1:maxit) {
      par = prox(t[i], theta - t[i] * df(theta))
      theta = par
      if(is.na(sqrt(sum((df(theta) ^ 2)))))
        stop("the subgradient vector is returning undefined value, change stepsize")
      if(sqrt(sum((df(theta) ^ 2))) < tol)
        break()
    }
  }
  value = fn(par)

  output = list("par" = as.vector(par), "value" = value, "iteration" = i)
  output
}

#' Coordinate descent optimisation for non smooth functions
#'
#' Function for minimising piecewise differentiable objective functions.
#'
#' @param theta Initial value.
#' @param fn Function to optimise.
#' @param theta_it Iterative theta for coordinate descent algorithm. Default
#' is NULL.
#' @param tol Relative convergence tolerance
#' @param maxit Maximum number of iterations if relative convergence is not
#' reached. Default is 1000.
#' @return The function returns a list with following components
#' \item{par}{The best set of parameters}
#' \item{value}{The value of the function corresponding to best set of parameters.}
#' \item{iteration}{Number of iterations taken}
#' @export

cdOptim = function(theta, fn, theta_it, tol = 1e-6, maxit = 1000) {
  for (j in 1:maxit) {
    theta.last = theta
    for(i in 1:length(theta)) {
      theta[i] = theta_it(theta, i)
    }
    if(sqrt(sum(((theta.last - theta) ^ 2))) < tol)
      break
  }
  value = fn(theta)

  output = list("par" = theta, "value" = value, "iteration" = j)
  output
}


#' Sequence of step sizes for optimization.
#'
#' Function to generate stepsize for sub-gradient optimization and proximal-gradient optimization
#' @param h Starting value. Default is 0.01
#' @param m Number of constant steps. Default is 1000
#' @param n Number of diminishing steps. Default is 1000
#' @return The sequence of stepsize
#' @export

stepsize = function(h = 0.01, m = 1000, n = 1000)
  if (n != 0) {
    c(rep(h, m), h / (1:n))
  } else {
    rep(h, m)
  }
