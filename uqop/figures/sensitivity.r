#install.packages("xtable")
#library(xtable)

# sensitivity analysis using weights

#gaia80a = read.table("stars.aps")
#gaia80b = read.table("stars.bprp")

#names(gaia80a) = c("ID", "Teff", "logg", "feh", "Av")
#names(gaia80b) = c("ID", paste("photon",1:86, sep=""))

#gaia80 = cbind(gaia80a, gaia80b[,2:87])

# with 86 predictors it's hard to explain from the plots

library(glmnet)
library(MCMCpack)
library(LPCM)

nsim = 1000

data(gaia)

ind = sample(1:nrow(gaia), 2000)

x = as.matrix(gaia[ind,5:20])
y1 = as.matrix(gaia[ind,4])

x = scale(x, scale = F)
y1 = scale(y1, scale = F)

wts  = 16 * rdirichlet(nsim, rep(1,16))

coef.gaia.wt = function(i, x, y, wt)
  gaia_coef = as.matrix(coef(cv.glmnet(x, y, penalty.factor = wt[i,], intercept = F)))

store = sapply(1:nsim, function(i) coef.gaia.wt(i, x, y1, wts))

boxplot(t(store[-1,]))
abline(h=0)

which(apply(store[-1,], 1, sd)==sort(apply(store[-1,], 1, sd))[16])

# cross-validation estimates
gaia_cv = cv.glmnet(x, y1, intercept = F)
gaia_coef = as.matrix(coef(gaia_cv)[-1])

plot(gaia_cv)

which(gaia_coef!=0)

gaia_cv2 = cv.glmnet(x[,-c(2)], y1, intercept = F)
gaia_coef2 = as.matrix(coef(gaia_cv2)[-1])

plot(gaia_cv2)

which(gaia_coef2!=0)

# correlation

correlation.gaia = cor(x,y1)
plot(1:16,correlation.gaia)
abline(h=0,lty=2)

# PCA
pca.gaia = princomp(x)

plot(pca.gaia)
