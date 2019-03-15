# R-Script for benchmarking UQOP Paris
#library(devtools)
#install_github("tathagatabasu/bootlasso")
#library(bootlasso)

nsim = 50

# system time change w.r.t no. of observations

sys_obs = function(n = 10){
  obs_sys_time_gen = function (i) {
    n = i * 10
    
    x = matrix(data = rnorm(3 * n), ncol = 3)
    b = as.matrix(c(-5, 0, 5))
    er = as.matrix(rnorm(n))
    y = x %*% b + er
    
    x = scale(x, scale = F)
    y = scale(y, scale = F)
    
    ts = opt_ts(0.005, 1000, 1000)
    wt = NULL
    
    system_time_sg_n = system.time(lasso_optim_sg(lambda = 1, x, y, ts, wt))[3]
    system_time_pg_n = system.time(lasso_optim_pg(lambda = 1, x, y, ts, wt))[3]
    system_time_cd_n = system.time(lasso_optim_cd(lambda = 1, x, y, 100, wt))[3]
    
    output = c("sg" = system_time_sg_n, "pg" = system_time_pg_n, "cd" = system_time_cd_n)
    
  }
  
  sys_time_obs = lapply(1:nsim, function(j) sapply(1:n, obs_sys_time_gen))
  
  system_time_sg_n = sapply(1:nsim, function(j) sys_time_obs[[j]][1,])
  system_time_pg_n = sapply(1:nsim, function(j) sys_time_obs[[j]][2,])
  system_time_cd_n = sapply(1:nsim, function(j) sys_time_obs[[j]][3,])
  
  system_time_sg_n_mean = rowMeans(system_time_sg_n)
  system_time_sg_n_sd = apply(system_time_sg_n, 1, sd)
  
  system_time_pg_n_mean = rowMeans(system_time_pg_n)
  system_time_pg_n_sd = apply(system_time_pg_n, 1, sd)
  
  system_time_cd_n_mean = rowMeans(system_time_cd_n)
  system_time_cd_n_sd = apply(system_time_cd_n, 1, sd)
  
  nobs = (1:n)*10
  
  yup_n = max(c(system_time_cd_n,system_time_pg_n,system_time_sg_n))
  
  plot(nobs, system_time_cd_n_mean, type = "l", ylim = c(0,yup_n), 
       ylab = "system time (sec)", xlab = "no. of observations",
       main = expression(paste("For, ", lambda, " = 1")), col = rgb(0,0,1))
  polygon(c(nobs, rev(nobs)), c((system_time_cd_n_mean-system_time_cd_n_sd),
                                      rev(system_time_cd_n_mean+system_time_cd_n_sd)),
          col = rgb(0,0,1,0.3), border = NA)
  
  lines(nobs, (system_time_pg_n_mean), col = 2)
  polygon(c(nobs, rev(nobs)), c((system_time_pg_n_mean-system_time_pg_n_sd),
                                      rev(system_time_pg_n_mean+system_time_pg_n_sd)),
          col = rgb(1,0,0,0.3), border = NA)
  
  lines(nobs, (system_time_sg_n_mean), col = 3)
  polygon(c(nobs, rev(nobs)), c((system_time_sg_n_mean-system_time_sg_n_sd),
                                      rev(system_time_sg_n_mean+system_time_sg_n_sd)),
          col = rgb(0,1,0,0.3), border = NA)
  
  output = list("sg" = system_time_sg_n, "pg" = system_time_pg_n, "cd" = system_time_cd_n)
  return(output)
}

sys_obs(200)
legend("topleft", legend=c("cd","pg", "sg"),
       col= c(4,2,3), lty=rep(1,3), cex=0.8)

# system time change w.r.t no. of predictors

sys_pred = function(p = 2){
  pred_sys_time_gen = function(i) {
    p = i + 2
    
    x = matrix(data = rnorm(100 * p), ncol = p)
    b = as.matrix(seq(-5, 5, length.out = p))
    er = as.matrix(rnorm(100))
    y = x %*% b + er
    
    x = scale(x, scale = F)
    y = scale(y, scale = F)
    
    ts = opt_ts(0.005, 1000, 1000)
    wt = NULL
    
    system_time_sg_p = system.time(lasso_optim_sg(lambda = 1, x, y, ts, wt))[3]
    system_time_pg_p = system.time(lasso_optim_pg(lambda = 1, x, y, ts, wt))[3]
    system_time_cd_p = system.time(lasso_optim_cd(lambda = 1, x, y, 100, wt))[3]
    
    output = c("sg" = system_time_sg_p, "pg" = system_time_pg_p, "cd" = system_time_cd_p)
    
  }
  
  sys_time_pred = lapply(1:nsim, function(j) sapply(1:p, pred_sys_time_gen))
  
  system_time_sg_p = sapply(1:nsim, function(j) sys_time_pred[[j]][1,])
  system_time_pg_p = sapply(1:nsim, function(j) sys_time_pred[[j]][2,])
  system_time_cd_p = sapply(1:nsim, function(j) sys_time_pred[[j]][3,])
  
  system_time_sg_p_mean = rowMeans(system_time_sg_p)
  system_time_sg_p_sd = apply(system_time_sg_p, 1, sd)
  
  system_time_pg_p_mean = rowMeans(system_time_pg_p)
  system_time_pg_p_sd = apply(system_time_pg_p, 1, sd)
  
  system_time_cd_p_mean = rowMeans(system_time_cd_p)
  system_time_cd_p_sd = apply(system_time_cd_p, 1, sd)
  
  npred = (1:p) + 2
  
  yup_p = max(c(system_time_cd_p,system_time_pg_p,system_time_sg_p))
  
  plot(npred, system_time_cd_p_mean, type = "l", ylim = c(0,yup_p), 
       ylab = "system time (sec)", xlab = "no. of predictors",
       main = expression(paste("For, ", lambda, " = 1")), col = rgb(0,0,1))
  polygon(c(npred, rev(npred)), c((system_time_cd_p_mean-system_time_cd_p_sd),
                                      rev(system_time_cd_p_mean+system_time_cd_p_sd)),
          col = rgb(0,0,1,0.3), border = NA)
  
  lines(npred, (system_time_pg_p_mean), col = 2)
  polygon(c(npred, rev(npred)), c((system_time_pg_p_mean-system_time_pg_p_sd),
                                      rev(system_time_pg_p_mean+system_time_pg_p_sd)),
          col = rgb(1,0,0,0.3), border = NA)
  
  lines(npred, (system_time_sg_p_mean), col = 3)
  polygon(c(npred, rev(npred)), c((system_time_sg_p_mean-system_time_sg_p_sd),
                                      rev(system_time_sg_p_mean+system_time_sg_p_sd)),
          col = rgb(0,1,0,0.3), border = NA)
  
  output = list("sg" = system_time_sg_p, "pg" = system_time_pg_p, "cd" = system_time_cd_p)
  return(output)
  
}

sys_pred(200)
#legend("topleft", legend=c("cd","pg", "sg"),
#       col= c(4,2,3), lty=rep(1,3), cex=0.8)

# total system time for solution path

sys_lambda = function(p = 3) {
  x = matrix(data = rnorm(100 * p), ncol = p)
  b = as.matrix(seq(-5, 5, length.out = p))
  er = as.matrix(rnorm(100))
  y = x %*% b + er
  
  x = scale(x, scale = F)
  y = scale(y, scale = F)
  
  lmax = max(abs(t(x) %*% y/diag(t(x) %*% x)))
  lambdas = as.matrix(exp(seq(-5, log(lmax), length.out = 101)))
  
  lambda_sys_time_gen = function(l) {
    ts = opt_ts(0.005, 1000, 1000)
    wt = NULL
    
    system_time_sg_l = system.time(lasso_optim_sg(lambda = l, x, y, ts, wt))[3]
    system_time_pg_l = system.time(lasso_optim_pg(lambda = l, x, y, ts, wt))[3]
    system_time_cd_l = system.time(lasso_optim_cd(lambda = l, x, y, 100, wt))[3]
    
    output = c("sg" = system_time_sg_l, "pg" = system_time_pg_l, "cd" = system_time_cd_l)
  }
  
  sys_time_lambda = lapply(1:nsim, function(j) sapply(lambdas, lambda_sys_time_gen))
  
  system_time_sg_l = sapply(1:nsim, function(j) sys_time_lambda[[j]][1,])
  system_time_pg_l = sapply(1:nsim, function(j) sys_time_lambda[[j]][2,])
  system_time_cd_l = sapply(1:nsim, function(j) sys_time_lambda[[j]][3,])
  
  system_time_sg_l_mean = rowMeans(system_time_sg_l)
  system_time_sg_l_sd = apply(system_time_sg_l, 1, sd)
  
  system_time_pg_l_mean = rowMeans(system_time_pg_l)
  system_time_pg_l_sd = apply(system_time_pg_l, 1, sd)
  
  system_time_cd_l_mean = rowMeans(system_time_cd_l)
  system_time_cd_l_sd = apply(system_time_cd_l, 1, sd)
  
  yup_l = max(c((system_time_cd_l_mean+system_time_cd_l_sd),
                (system_time_pg_l_mean+system_time_pg_l_sd),
                (system_time_sg_l_mean+system_time_sg_l_sd)))
  
  plot(log(lambdas), (system_time_cd_l_mean), type = "l", ylim = c(0,yup_l), 
       ylab = "system time (sec)", xlab = expression(paste("log",lambda)),
       main = paste("no of predictors = ", p), col = rgb(0,0,1))
  polygon(c(log(lambdas), rev(log(lambdas))), c((system_time_cd_l_mean-system_time_cd_l_sd),
                                      rev(system_time_cd_l_mean+system_time_cd_l_sd)),
          col = rgb(0,0,1,0.3), border = NA)
  
  lines(log(lambdas), (system_time_pg_l_mean), col = 2)
  polygon(c(log(lambdas), rev(log(lambdas))), c((system_time_pg_l_mean-system_time_pg_l_sd),
                                      rev(system_time_pg_l_mean+system_time_pg_l_sd)),
          col = rgb(1,0,0,0.3), border = NA)
  
  lines(log(lambdas), (system_time_sg_l_mean), col = 3)
  polygon(c(log(lambdas), rev(log(lambdas))), c((system_time_sg_l_mean-system_time_sg_l_sd),
                     rev(system_time_sg_l_mean+system_time_sg_l_sd)),
          col = rgb(0,1,0,0.3), border = NA)
  
  output = list("sg" = system_time_sg_l, "pg" = system_time_pg_l, "cd" = system_time_cd_l)
  return(output)
}

sys_lambda(3)
#legend("topright",bg="transparent", legend=c("cd","pg", "sg"),
#       col= c(4,2,3), lty=rep(1,3), cex=0.8)
sys_lambda(9)
#legend("right", legend=c("cd","pg", "sg"),
#       col= c(4,2,3), lty=rep(1,3), cex=0.8)
sys_lambda(27)
#legend("topright", legend=c("cd","pg", "sg"),
#       col= c(4,2,3), lty=rep(1,3), cex=0.8)
sys_lambda(81)
#legend("topright", legend=c("cd","pg", "sg"),
#       col= c(4,2,3), lty=rep(1,3), cex=0.8)


# convergence rate

conv_track = function(p = 3){
  n_it = 2000
  ts = opt_ts(0.005, n_it, 0)
  wt = NULL
  
  # function over-ride
  
  lasso_optim_cd2 = function (lambda, x, y, n_it = 100, wt) 
  {
    if ((is.null(wt) == T) | (length(wt) != ncol(x))) 
      wt = rep(1, ncol(x))
    else wt = ncol(x) * wt/sum(wt)
    beta0 = as.matrix(rep(0, ncol(x)))
    s = soft(lambda, wt)
    f = function(beta) square_lasso_f(lambda, x, y, beta, wt)
    v = function(i, beta) st_f(i, x, y, beta)
    x_it = function(beta, i) s(v(i, beta), i)
    cd_optim2(x = beta0, f = f, x_it = x_it, n_it = n_it)
  }
  
  cd_optim2=function (x, f, x_it, n_it = 100) 
  {
    track = c()
    for (j in 1:n_it) {
      x.last = x
      for (i in 1:length(x)) {
        x[i] = x_it(x, i)
      }
      track = cbind(track, sum(abs(x.last - x)))
      if (sum(abs(x.last - x)) < 1e-07) 
        break
    }
    track
  }
  
  lasso_optim_pg2 = function (lambda, x, y, ts, wt) 
  {
    if ((is.null(wt) == T) | (length(wt) != ncol(x))) 
      wt = rep(1, ncol(x))
    else wt = ncol(x) * wt/sum(wt)
    beta0 = as.matrix(rep(0, ncol(x)))
    f = function(beta) square_f(x, y, beta)
    df = function(beta) square_df(x, y, beta)
    pg = lasso_p(lambda, wt)
    pg_optim2(x = beta0, f = f, df = df, pg = pg, ts = ts)
  }
  
  pg_optim2 = function (x, f, df, pg, ts) 
  {
    track = c()
    for (t in ts) {
      x.last = x
      x = pg(t, x - t * df(x))
      track = cbind(track, sum(abs(x.last - x)))
    }
    track
  }
  
  lasso_optim_sg2 = function (lambda, x, y, ts, wt) 
  {
    if ((is.null(wt) == T) | (length(wt) != ncol(x))) 
      wt = rep(1, ncol(x))
    else wt = ncol(x) * wt/sum(wt)
    beta0 = as.matrix(rep(0, ncol(x)))
    f = function(beta) square_lasso_f(lambda, x, y, beta, wt)
    df = function(beta) square_lasso_df(lambda, x, y, beta, wt)
    sg_optim2(x = beta0, f = f, df = df, ts = ts)
  }
  
  sg_optim2 = function (x, f, df, ts) 
  {
    fx = f(x)
    x.best = x
    track = c()
    fx.best = fx
    for (t in ts) {
      x.last = x
      x = x - t * df(x)
      fx = f(x)
      if (fx < fx.best) {
        x.best = x
        fx.best = fx
      }
      track = cbind(track, sum(abs(x.last - x)))
    }
    track
  }
  
  
  sg_sol_track = function(j){
    x = matrix(data = rnorm(100 * p), ncol = p)
    b = as.matrix(seq(-.5, .5, length.out = p))
    er = as.matrix(rnorm(100))
    y = x %*% b + er
    
    x = scale(x, scale = F)
    y = scale(y, scale = F)
    lm_coef = lm(y~x-1)$coef
    
    sol_sg = lasso_optim_sg2(lambda = 1, x, y, ts, wt)
    
    return(sol_sg)
  }
  
  pg_sol_track = function(j){
    x = matrix(data = rnorm(100 * p), ncol = p)
    b = as.matrix(seq(-.5, .5, length.out = p))
    er = as.matrix(rnorm(100))
    y = x %*% b + er
    
    x = scale(x, scale = F)
    y = scale(y, scale = F)
    lm_coef = lm(y~x-1)$coef
    
    sol_pg = lasso_optim_pg2(lambda = 1, x, y, ts, wt)
    
    return(sol_pg)
  }
  
  cd_sol_track = function(j){
    x = matrix(data = rnorm(100 * p), ncol = p)
    b = as.matrix(seq(-.5, .5, length.out = p))
    er = as.matrix(rnorm(100))
    y = x %*% b + er
    
    x = scale(x, scale = F)
    y = scale(y, scale = F)
    lm_coef = lm(y~x-1)$coef
    
    sol_cd = lasso_optim_cd2(lambda = 1, x, y, n_it = n_it, wt)
    sol_cd[(length(sol_cd)+1):n_it] = 0
    
    return(sol_cd)
  }
  
  
  conv_track_cd = sapply(1:nsim, cd_sol_track)
  conv_track_sg = sapply(1:nsim, sg_sol_track)
  conv_track_pg = sapply(1:nsim, pg_sol_track)
  
  conv_track_cd_mean = rowMeans(conv_track_cd)
  conv_track_cd_sd = apply(conv_track_cd, 1, sd)
  
  conv_track_sg_mean = rowMeans(conv_track_sg)
  conv_track_sg_sd = apply(conv_track_sg, 1, sd)
  
  conv_track_pg_mean = rowMeans(conv_track_pg)
  conv_track_pg_sd = apply(conv_track_pg, 1, sd)
  
  yup_it = max(c((conv_track_cd_mean+conv_track_cd_sd),
                 (conv_track_pg_mean+conv_track_pg_sd),
                 (conv_track_sg_mean+conv_track_sg_sd)))
  
  plot(log10(1:n_it), conv_track_cd_mean, type = "l", ylim = c(0,yup_it), 
       ylab = "absolute error", xlab = "no. of iterations in power of 10",
       main = expression(paste("For, ", lambda, " = 1")), col = rgb(0,0,1))
  polygon(c(log10(1:n_it), rev(log10(1:n_it))), c((conv_track_cd_mean-conv_track_cd_sd),
                                              rev(conv_track_cd_mean+conv_track_cd_sd)),
          col = rgb(0,0,1,0.3), border = NA)
  
  lines(log10(1:n_it), (conv_track_pg_mean), col = 2)
  polygon(c(log10(1:n_it), rev(log10(1:n_it))), c((conv_track_pg_mean-conv_track_pg_sd),
                                              rev(conv_track_pg_mean+conv_track_pg_sd)),
          col = rgb(1,0,0,0.3), border = NA)
  
  lines(log10(1:n_it), (conv_track_sg_mean), col = 3)
  polygon(c(log10(1:n_it), rev(log10(1:n_it))), c((conv_track_sg_mean-conv_track_sg_sd),
                                              rev(conv_track_sg_mean+conv_track_sg_sd)),
          col = rgb(0,1,0,0.3), border = NA)
  
  # zoom 
  
  plot(log10(1:n_it), conv_track_cd_mean, type = "l", 
       ylim = c(0,yup_it / (max(conv_track_cd_mean) / 
                              max(c(conv_track_pg_mean, conv_track_sg_mean)))), 
       ylab = "absolute error", xlab = "no. of iterations in power of 10",
       main = expression(paste("For, ", lambda, " = 1 (zoomed)")), col = rgb(0,0,1))
  polygon(c(log10(1:n_it), rev(log10(1:n_it))), c((conv_track_cd_mean-conv_track_cd_sd),
                                                  rev(conv_track_cd_mean+conv_track_cd_sd)),
          col = rgb(0,0,1,0.3), border = NA)
  
  lines(log10(1:n_it), (conv_track_pg_mean), col = 2)
  polygon(c(log10(1:n_it), rev(log10(1:n_it))), c((conv_track_pg_mean-conv_track_pg_sd),
                                                  rev(conv_track_pg_mean+conv_track_pg_sd)),
          col = rgb(1,0,0,0.3), border = NA)
  
  lines(log10(1:n_it), (conv_track_sg_mean), col = 3)
  polygon(c(log10(1:n_it), rev(log10(1:n_it))), c((conv_track_sg_mean-conv_track_sg_sd),
                                                  rev(conv_track_sg_mean+conv_track_sg_sd)),
          col = rgb(0,1,0,0.3), border = NA)
  
  output = list("sg" = conv_track_sg, "pg" = conv_track_pg, "cd" = conv_track_cd)
  return(output)
}

#conv_track(3)
#legend("topright", legend=c("cd","pg", "sg"),
#       col= c(4,2,3), lty=rep(1,3), cex=0.8)

library(lattice)
library(gridExtra) # grid.arrange

contplot = function(a, b) {
  xs = seq(-1,1,0.01)
  grid = expand.grid(x=xs, y=xs)
  grid$z = (a * (abs(grid$x) ^ 1) + (b * abs(grid$y) ^ 1)) ^ (1/1)
  cols = gray(seq(0.2, 1, 0.01))
  contourplot(
    z~x*y, grid,
    labels=FALSE, xlab=NULL, ylab=NULL,
    aspect=1, scales=list(draw=FALSE))
}

pdf("images/lqplots.pdf", width=12, height=3)
grid.arrange(
  contplot(1,1),
  contplot(1,3),
  ncol=2)
dev.off()
