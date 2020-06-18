# spBayes example code
library(coda)
library(spBayes)

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