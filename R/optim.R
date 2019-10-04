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
#' @param method Method to be used for optimisation. Default is subgradient
#' method.
#' @param tol Relative convergence tolerance
#' @param maxit Maximum number of iterations if relative convergence is not
#' reached. Default is 1000.
#' @param h Stepsize for updating the set of parameters. Default is 0.001.
#' @param prox Proximal operator for proximal gradient method. Default is 
#' NULL.
#' @param acl Flag for accelerated proximal gradient method. Default is false.
#' @param theta_it Iterative theta for coordinate descent algorithm. Default
#' is NULL.
#' @return The function returns a list with following components
#' \item{par}{The best set of parameters}
#' \item{value}{The value of the function corresponding to best set of parameters.}
#' \item{iteration}{Number of iterations taken}
#' \item{method}{The method used for the optimisation}
#' @export

nonsmoothOptim = function(theta, fn, df, method = "SG", tol = 1e-5, maxit = 1000, h = 0.001,
                        prox = NULL, acl = FALSE, theta_it = NULL)
{
  if(length(theta) != length(df(theta)))
    stop("Number of parameters does not match with length of the sub gradient vector")
  
  if(is.null(prox) & any(method == c("PG", "pg"))){
    warning("No proximal operator provided, the method defaults to sub gradient method")
    method = "SG"
  }
    
  if(is.null(theta_it) && any(method == c("CD", "cd"))){
    warning("No iterative theta provided, the method defaults to sub gradient method")
    method = "SG"
  }
  
  # Sub-gradient optimisation.
  sg = function(theta, fn, df, tol, maxit, h) {
    value = fn(theta)
    par = theta.last = theta
    t = c(rep(h, maxit %/% 2), h / (1:(maxit - maxit %/% 2)))
    for (i in 1:maxit) {
      theta = theta.last - t[i] * df(theta.last)
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
    
    output = list("par" = par, "value" = value, "iteration" = i, "method" = "sub gradient")
    output
  }
  
  # Proximal gradient optimisation.
  
  pg = function(theta, fn, df, prox, acl, tol, maxit, h) {
    par = theta.prev = theta
    t = c(rep(h, maxit %/% 2), h / (1:(maxit - maxit %/% 2)))
    if(isTRUE(acl)){
      for (i in 1:maxit) {
        v = theta + (i - 2) * (theta - theta.prev) / (i + 1)
        par = prox(t[i], v - t[i] * df(v))
        theta.prev = theta
        theta = par
        if(sqrt(sum((df(par) ^ 2))) < tol)
          break()
      }
      throw = "accelerated proximal gradient"
    }
    else{
      for (i in 1:maxit) {
        par = prox(t[i], theta - t[i] * df(theta))
        theta = par
        if(sqrt(sum((df(par) ^ 2))) < tol)
          break()
      }
      throw = "proximal gradient"
    }
    value = fn(par)
    
    output = list("par" = par, "value" = value, "iteration" = i, "method" = throw)
    output
  }
  
  # Co-ordinate descent optimisation.
  
  cd = function(theta, fn, theta_it, tol, maxit) {
    for (j in 1:maxit) {
      theta.last = theta
      for(i in 1:length(theta)) {
        theta[i] = theta_it(theta, i)
      }
      if(sqrt(sum(((theta.last - theta) ^ 2))) < tol)
        break
    }
    value = fn(theta)
    
    output = list("par" = theta, "value" = value, "iteration" = j, "method" = "coordinate descent")
    output
  }
  
  if(any(method == c("SG", "sg"))) 
    output = sg(theta = theta, fn = fn, df = df, tol = tol, maxit = maxit, h = h)
  if(any(method == c("PG", "pg")))
    output = pg(theta = theta, fn = fn, df = df, prox = prox, acl = acl, 
                tol = tol, maxit = maxit, h = h)
  if(any(method == c("CD", "cd")))
    output = cd(theta = theta, fn = fn, theta_it = theta_it, tol = tol, maxit = maxit)
  
  
  output
}
