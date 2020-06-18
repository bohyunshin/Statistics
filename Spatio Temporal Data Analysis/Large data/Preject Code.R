rm(list=ls())

#####################################################################################
###################### import libraries and load data  ##############################
#####################################################################################
library(dplyr)
library(spNNGP)
library(spBayes)
library(fields)
library(LatticeKrig)
library(INLA)
library('sp')
library(FRK)
library(dggrids)
library(gridExtra)
library(ggplot2)
library(tidyverse)
load('/Users/shinbo/Desktop/20-1/spatial_statistics/high_dimension problems/comparing_paper/Data/AllSatelliteTemps.RData')
fig.dir = '/Users/shinbo/Desktop/20-1/spatial_statistics/high_dimension problems/tex/figure/'

#####################################################################################
########################## define function to evaluate ##############################
#####################################################################################
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



#####################################################################################
############################### 1. prepare data  ####################################
#####################################################################################

# na observation in sat data
sum(is.na(all.sat.temps$TrueTemp))

# exclude missing data
all.sat.temps = all.sat.temps %>% filter(! is.na(TrueTemp))

# split satellite data
sat.train = all.sat.temps %>% filter(! is.na(MaskTemp) & ! is.na(TrueTemp))
sat.test = all.sat.temps %>% filter( is.na(MaskTemp)  & ! is.na(TrueTemp))
dim(sat.train)[1] ; dim(sat.test)[1]

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

#####################################################################################
################################ 2. SNNGP Modeling  #################################
#####################################################################################
n.samples <- 5000
theta.alpha <- as.matrix(expand.grid(seq(0.1,1,length.out=15), seq(3/10,3/0.1,length.out=15)))
colnames(theta.alpha) <- c("alpha","phi")

start = Sys.time()
result.snngp <- spConjNNGP(TrueTemp ~ 1, coords=c('Lon','Lat'), data=sat.train, cov.model="exponential", sigma.sq.IG=c(2, 40),
                           n.neighbors=10, theta.alpha=theta.alpha,
                           k.fold = 2, score.rule = "crps", fit.rep=TRUE, n.samples=n.samples, n.omp.threads=18,
                           verbose=TRUE)
modeling.time.snngp = start - Sys.time()

head(theta.alpha)
result.snngp$theta.alpha
head(result.snngp$k.fold.scores)
summary(result.snngp)

############ krigging ############
x.val = cbind(rep(1, nrow(sat.test)) )
y.act = sat.test$TrueTemp
sat.coords = as.matrix(sat.test[,1:2])

start = Sys.time()
sat.pred = predict(result.snngp, X.0 = x.val, coords.0 = sat.coords, n.omp.threads=1)
krig.time.snngp = Sys.time() - start

sat.test$snngp_fitted = sat.pred$y.0.hat
sat.test$snngp_fitted_sd = sqrt(sat.pred$y.0.hat.var)

############ plot result ############
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

############ evaluation ############
evaluate(sat.test$TrueTemp, sat.test$snngp_fitted, sat.test$snngp_fitted_sd)

############ computation time  ############
modeling.time.snngp + krig.time.snngp



#####################################################################################
################################# 2. INLA Modeling  #################################
#####################################################################################

############ define mesh  ############
# mesh must cover the entire spatial domain of interest
mesh = inla.mesh.2d(loc.domain = as.matrix(sat.train[,1:2]), cutoff = 0.05, offset = c(0.1, 0.4),  max.edge = c(0.05, 0.5))

############ define spde model  ############
# just matern or pcmatern?
# when done with pcmatern, the results were bad, so go with the ordinary matern
spde.matern = inla.spde2.matern(mesh = mesh, alpha = 2) # alpha = 2 is default value
spde.pcmatern = inla.spde2.pcmatern(mesh = mesh, alpha = 2,
                                    prior.range = c(0.3,0.5),
                                    prior.sigma = c(10, 0.01))

############ define projection matrix  ############
# The basis matrix A in INLA paper
coords.train = as.matrix(sat.train[,1:2])
A.train = inla.spde.make.A(mesh, loc = coords.train)
# mesh$n: 4059 basis
#dim(A.train): spatial process of each observation consists of 4059 basis

############ data stack  ############
stack.train = inla.stack(
  data = list(y = sat.train$True),
  A = list(A.train, 1),
  effects = list(i = 1:spde.matern$n.spde,
                 beta0 = rep(1, nrow(sat.train))),
  tag = 'sat.train'
)

start = Sys.time()
result = inla(y ~ 0 + beta0 + f(i, model = spde.matern),
              data = inla.stack.data(stack.train),
              control.predictor = list(A = inla.stack.A(stack.train)), control.compute=list(config = TRUE))
modeling.time.inla = Sys.time() - start

# summary of beta0
result$summary.fixed

# summary of the other parameters
result$summary.hyperpar

# marginal plot for beta0
plot(result$marginals.fixed[[1]], ty='l', xlab = expression(beta[0]), ylab = 'Density')
result$summary.fitted.values



############ krigging  ############

############ stack test data ############
coords.test = as.matrix(sat.test[,1:2])
A.test = inla.spde.make.A(mesh, loc = coords.test)
stack.test = inla.stack(
  data = list(y = NA),
  A = list(A.test, 1),
  effects = list(i = 1:spde.matern$n.spde,
                 beta0 = rep(1, nrow(sat.test))),
  tag = 'sat.test'
)

############ full stack ############
stack.full = inla.stack(stack.train, stack.test)
start = Sys.time()
pred = inla(y ~ 0 + beta0 + f(i, model = spde.matern),
            data = inla.stack.data(stack.full),
            control.mode = list(theta = result$mode$theta, restart = FALSE),
            control.predictor = list(A = inla.stack.A(stack.full), compute=TRUE)
)
krig.time.inla = Sys.time() - start

############ krigging results ############
index.pred = inla.stack.index(stack.full, 'sat.test')$data
krig.mean = pred$summary.fitted.values[index.pred, 'mean']
krig.sd = pred$summary.fitted.values[index.pred, 'sd']

sat.test$inla_fitted = krig.mean
sat.test$inla_fitted_sd = krig.sd

############ plot result ############
ggplot(sat.test, aes(Lon,Lat)) + 
  theme_bw() + 
  geom_raster(aes(fill = inla_fitted_sd)) +
  scale_fill_distiller(palette="Spectral", name="", limits = c(0,max_sd))
ggsave("/Users/shinbo/Desktop/20-1/spatial_statistics/high_dimension problems/tex/figure/INLA krigging sd.png")

ggplot(sat.test, aes(Lon,Lat)) + 
  theme_bw() + 
  geom_raster(aes(fill = inla_fitted)) +
  scale_fill_distiller(palette="Spectral", name="", limits = c(min_krig, max_krig))
#+ ggtitle("Krigging result") + 
#theme(plot.title = element_text(hjust = 0.5))
ggsave("/Users/shinbo/Desktop/20-1/spatial_statistics/high_dimension problems/tex/figure/INLA krigging.png")

############ evaluate performance ############
evaluate(sat.test$TrueTemp, sat.test$inla_fitted, sat.test$inla_fitted_sd)

############ evaluate computation time ############
modeling.time.inla + krig.time.inla

#####################################################################################
###################################### 4. FRK  ######################################
#####################################################################################

tol = 0.01
f = TrueTemp ~ 1

############ manually assign BAUs ############
BAUs <- all.sat.temps                    
BAUs$Missing <- is.na(BAUs$MaskTemp)      
BAUs$MaskTemp <- NULL                     
BAUs$TrueTemp <- NULL
BAUs$fs <- 1
# convert to SpatialPointsDataFrame
coordinates(BAUs)  <- ~Lon+Lat        
gridded(BAUs) <- TRUE  


train.FRK <- subset(all.sat.temps,!is.na(MaskTemp))     
coordinates(train.FRK)  <- ~Lon+Lat

############ create basis ############
# hyperparameters were chosen according to 
# paper A Case Study Competition Among Methods for Analyzing Large Spatial Data
basis <- auto_basis(plane(),          
                    data = train.FRK,       
                    regular = 0,      
                    nres = 3,         
                    scale_aperture = 1) 
start = Sys.time()
S <- FRK(f = f,                       
         data = train.FRK,                  
         basis = basis,               
         BAUs = BAUs,                 
         tol = tol)
modeling.time.frk = Sys.time() - start

############ krig ############
start = Sys.time()
BAUs_pred <- predict(S)           # predict over all BAUs
krig.time.frk = Sys.time() - start

BAUs_pred_df <- data.frame(BAUs_pred) # convert to data frame

# select the test data
tmp = BAUs_pred_df %>% filter(Missing == T)

sat.test$frk_fitted = tmp$mu
sat.test$frk_fitted_sd = tmp$sd

############ plot result ############
ggplot(sat.test, aes(Lon,Lat)) + 
  theme_bw() + 
  geom_raster(aes(fill = frk_fitted)) +
  scale_fill_distiller(palette="Spectral", name="", limits = c(min_krig, max_krig))
#+ ggtitle("Plot for Pred Y")
ggsave(paste(fig.dir, "frk krigging.png",sep=''))

ggplot(sat.test, aes(Lon,Lat)) + 
  theme_bw() + 
  geom_raster(aes(fill = frk_fitted_sd)) +
  scale_fill_distiller(palette="Spectral", name="", limits = c(0, max_sd))
#+ ggtitle("Plot for Pred Y")
ggsave(paste(fig.dir, "frk krigging sd.png",sep=''))

############ evaluate performance ############
evaluate(sat.test$TrueTemp, sat.test$frk_fitted, sat.test$frk_fitted_sd)

############ evaluate computation time ############
modeling.time.frk + krig.time.frk

#####################################################################################
########################### 5. Predictive Process Model  ############################
#####################################################################################

#####################################################################################
# Note: The Model was not completely finished due to the excessive computation time
# So I just attached the related code and do not report the results
#####################################################################################

knots = as.matrix(read.table('/Users/shinbo/Desktop/20-1/spatial_statistics/high_dimension problems/knots'))
names(knots) = c('Lon', 'Lat')
p = 1
n.samples = 5000

############ setting priors ############
starting <- list("phi"=3/0.5, "sigma.sq"=50, "tau.sq"=1) 
tuning <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1)
priors.1 <- list("beta.Norm"=list(rep(0,p), diag(1000,p)), "phi.Unif"=c(3/1, 3/0.1), "sigma.sq.IG"=c(2, 2),
                 "tau.sq.IG"=c(2, 0.1))
coords = as.matrix(sat.train[,c('Lon','Lat')])

ppm = spLM(TrueTemp ~ 1, data = sat.train, coords = coords, knots = knots, starting = starting,
           tuning = tuning, priors = priors.1, cov.model = cov.model, modified.pp = T,
           n.samples = n.samples)

coords.test = as.matrix(sat.test[,c('Lon','Lat')])

# this krigging part takes more than one day.
# due to the time limit, I just decided not to include the results of ppm
ppm.pred = spPredict(ppm, pred.covars = as.matrix(rep(1, dim(sat.test)[1])), 
                     pred.coords = coords.test)
# if krigging is done appropriately,
y.hat <- apply(m.1.pred$p.y.predictive.samples, 1, quant)

save.image('/Users/shinbo/Desktop/20-1/spatial_statistics/high_dimension problems/final_workspace.RData')
