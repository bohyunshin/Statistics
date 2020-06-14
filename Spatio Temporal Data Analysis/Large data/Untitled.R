load('/Users/shinbo/Desktop/20-1/spatial_statistics/high_dimension problems/comparing_paper/Data/AllSimulatedTemps.RData')



library(spNNGP)
library(dplyr)
library(MBA)
library(fields)

data('BCEF')
BCEF.mod = BCEF %>% filter(holdout == 0)
n.samples <- 5000
starting <- list("phi"=3/2, "sigma.sq"=40, "tau.sq"=1)
priors <- list("phi.Unif"=c(3/10, 3/0.1), "sigma.sq.IG"=c(2, 40), "tau.sq.IG"=c(2, 10))
cov.model <- "exponential"
tuning <- list("phi"=0.02)
bcef.s <- spNNGP(FCH~PTC, coords=c("x","y"), data=BCEF.mod, starting=starting, method="latent", n.neighbors=10,
                 tuning=tuning, priors=priors, cov.model=cov.model, n.samples=n.samples, n.omp.threads=18, n.report=2500, fit.rep=TRUE, sub.sample=list(start=4000, thin=10))


tuning <- list("phi"=0.01, "sigma.sq"=0.01, "tau.sq"=0.005)
bcef.r <- spNNGP(FCH~PTC, coords=c("x","y"), data=BCEF.mod, starting=starting, method="response", n.neighbors=10,
                 tuning=tuning, priors=priors, cov.model=cov.model, n.samples=n.samples, 
                 n.omp.threads=18, n.report=2500, fit.rep=TRUE, sub.sample=list(start=4000, thin=10), verbose=FALSE)

theta.alpha <- as.matrix(expand.grid(seq(0.1,1,length.out=15), seq(3/10,3/0.1,length.out=15)))
colnames(theta.alpha) <- c("alpha","phi")
bcef.c <- spConjNNGP(FCH~PTC, coords=c("x","y"), data=BCEF.mod, cov.model="exponential", sigma.sq.IG=c(2, 40),
                     n.neighbors=10, theta.alpha=theta.alpha,
                     k.fold = 2, score.rule = "crps", fit.rep=TRUE, n.samples=200, n.omp.threads=18,
                     verbose=FALSE)

summary(bcef.c)
spDiag(bcef.c)

crps.surf <- mba.surf(bcef.c$k.fold.scores[,c("phi","alpha","crps")], no.X=100, no.Y=100)$xyz.est

image.plot(crps.surf, xlab="phi", ylab="alpha=tau^2/sigma^2", main="CRPS (lower is better)")

# prediction
theta.alpha <- as.vector(bcef.c$theta.alpha)
names(theta.alpha) <- c("phi", "alpha")
