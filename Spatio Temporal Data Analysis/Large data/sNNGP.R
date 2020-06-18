

# basic libraries to handle data
library(dplyr)
# libraries to model spatial data
library(spNNGP)
library(spBayes)
library(fields)
library(LatticeKrig)

# for visualization
library(gridExtra)
library(ggplot2)
library(tidyverse)

# load data
load('/Users/shinbo/Desktop/20-1/spatial_statistics/high_dimension problems/comparing_paper/Data/AllSatelliteTemps.RData')
load('/Users/shinbo/Desktop/20-1/spatial_statistics/high_dimension problems/comparing_paper/Data/AllSimulatedTemps.RData')
fig.dir = '/Users/shinbo/Desktop/20-1/spatial_statistics/high_dimension problems/tex/figure/'

rmse = function(act, pred){
  return(sqrt(mean((act - pred)^2)))
}
mae = function(act, pred){
  return(mean(abs(act - pred)))
}
mycrps <- function(pred, pred.se, act) {
  z <- as.numeric((act - pred) / pred.se)
  scores <- pred.se * (z *(2 * pnorm(z, 0, 1) - 1) +
                         2 * dnorm(z, 0, 1) - 1/sqrt(pi))
  return(mean(scores))
}

evaluate = function(act, pred, pred.se){
  cat('rmse, mae, crps scores')
  return( c(rmse(act,pred), mae(act,pred), mycrps(pred,pred.se,act)) )
}

## int, cvg는 아직.

############ prepare data ############

# na observation in sat data
sum(is.na(all.sat.temps$TrueTemp))

# because simulated data, no na observation in sim.data
sum(is.na(all.sim.data$TrueTemp))

# exclude missing data
all.sat.temps = all.sat.temps %>% filter(! is.na(TrueTemp))

# split satellite data
sat.train = all.sat.temps %>% filter(! is.na(MaskTemp) & ! is.na(TrueTemp))
sat.test = all.sat.temps %>% filter( is.na(MaskTemp)  & ! is.na(TrueTemp))
dim(sat.train)[1] ; dim(sat.test)[1]

# split sim data
sim.train = all.sim.data %>% filter(! is.na(MaskTemp))
sim.test = all.sim.data %>% filter( is.na(MaskTemp))
dim(sim.train)[1] ; dim(sim.test)[1]

knots = read.table('/Users/shinbo/Desktop/20-1/spatial_statistics/high_dimension problems/knots')

# plot for all data
ggplot(all.sat.temps, aes(Lon,Lat)) + 
  theme_bw() + 
  geom_raster(aes(fill = TrueTemp)) +
  scale_fill_distiller(palette="Spectral", name="") + 
  geom_point(data=knots, aes(V1, V2))
  #+ ggtitle("Plot for observed temperature")
ggsave(paste(fig.dir, 'temperature all.png', sep = ''))

# plot test data
ggplot(sat.test, aes(Lon,Lat)) + 
  theme_bw() + 
  geom_raster(aes(fill = TrueTemp)) +
  scale_fill_distiller(palette="Spectral", name="")
#+ ggtitle("Plot for True Y")
ggsave(paste(fig.dir, 'temperature test.png', sep = ''))

############ sNNGP ############
n.samples <- 5000
theta.alpha <- as.matrix(expand.grid(seq(0.1,1,length.out=15), seq(3/10,3/0.1,length.out=15)))
colnames(theta.alpha) <- c("alpha","phi")
result.snngp <- spConjNNGP(TrueTemp ~ 1, coords=c('Lon','Lat'), data=sat.train, cov.model="exponential", sigma.sq.IG=c(2, 40),
                     n.neighbors=10, theta.alpha=theta.alpha,
                     k.fold = 2, score.rule = "crps", fit.rep=TRUE, n.samples=n.samples, n.omp.threads=18,
                     verbose=TRUE)

head(theta.alpha)
result.snngp$theta.alpha
head(result.snngp$k.fold.scores)
summary(result.snngp)

# run time
result.snngp$run.time

# krigging
x.val = cbind(rep(1, nrow(sat.test)) )
y.act = sat.test$TrueTemp
sat.coords = as.matrix(sat.test[,1:2])
sat.pred = predict(result.snngp, X.0 = x.val, coords.0 = sat.coords, n.omp.threads=1)
sat.pred$run.time
sat.test$snngp_fitted = sat.pred$y.0.hat
sat.test$snngp_fitted_sd = sqrt(sat.pred$y.0.hat.var)

############ plot result ############

# plot for true yhat: predicted little lower?
max_krig = max(sat.test$snngp_fitted, sat.test$inla_fitted, sat.test$frk_fitted)
min_krig = min(sat.test$snngp_fitted, sat.test$inla_fitted, sat.test$frk_fitted)
max_sd = max(sat.test$snngp_fitted_sd, sat.test$inla_fitted_sd, sat.test$frk_fitted_sd)

ggplot(sat.test, aes(Lon,Lat)) + 
  theme_bw() + 
  geom_raster(aes(fill = snngp_fitted)) +
  scale_fill_distiller(palette="Spectral", name="", limits = c(min_krig, max_krig))
  #+ ggtitle("Plot for Pred Y")
ggsave(paste(fig.dir, "snngp krigging.png",sep=''))

ggplot(sat.test, aes(Lon,Lat)) + 
  theme_bw() + 
  geom_raster(aes(fill = snngp_fitted_sd)) +
  scale_fill_distiller(palette="Spectral", name="", limits = c(0, max_sd))
#+ ggtitle("Plot for Pred Y")
ggsave(paste(fig.dir, "snngp krigging sd.png",sep=''))

#  evaluation
evaluate(sat.test$TrueTemp, sat.test$snngp_fitted, sat.test$snngp_fitted_sd)
