library(INLA)
m = 50
points = matrix(runif(m*2), m, 2)
mesh <- inla.mesh.2d(loc = points, cutoff = 0.05, offset = c(0.1, 0.4), max.edge = c(0.05, 0.5))
plot(mesh)
bnd <- inla.nonconvex.hull(points, convex = 0.12)
mesh <- inla.mesh.2d(boundary = bnd, cutoff = 0.05, max.edge = c(0.1))
plot(mesh)

inla.spde.models()
spde <- inla.spde2.matern(mesh, alpha = 2)

m = 100
points = matrix(runif(m*2),m,2)
mesh = inla.mesh.create.helper(
  points=points, cutoff=0.05, offset=c(0.1,0.4), max.edge=c(0.05,0.5) )
plot(mesh)
points(points[,1],points[,2])


A=inla.spde.make.A( mesh,
                    loc=points,
                    index=rep(1:m,times=2),
                    repl=rep(1:2,each=m) )
A[3,]

library(gstat)
library(INLA)
library(maptools)
data('meuse')
coordinates(meuse) <- ~x+y
proj4string(meuse) <- CRS("+init=epsg:28992")
gridded(meuse) = TRUE

#Code from gstat
data(meuse.grid)
coordinates(meuse.grid) = ~x+y
proj4string(meuse.grid) <- CRS("+init=epsg:28992")
gridded(meuse.grid) = TRUE
dim(meuse)

#Boundary
meuse.bdy <- unionSpatialPolygons(
  as(meuse.grid, "SpatialPolygons"), rep (1, length(meuse.grid))
)
meuse.grid$
#Define mesh
pts <- meuse.bdy@polygons[[1]]@Polygons[[1]]@coords
mesh <- inla.mesh.2d(loc.domain = pts, max.edge = c(150, 500),
                     offset = c(100, 250) )
par(mar = c(0, 0, 0, 0))
plot(mesh, asp = 1, main = "")
lines(pts, col = 3, with = 2)

#Create SPDE
meuse.spde <- inla.spde2.matern(mesh = mesh, alpha = 2)
A.meuse <- inla.spde.make.A(mesh = mesh, loc = coordinates(meuse))
s.index <- inla.spde.make.index(name = "spatial.field",
                                n.spde = meuse.spde$n.spde)


#Create data structure
meuse.stack <- inla.stack(data  = list(zinc = meuse$zinc),
                          A = list(A.meuse, 1),
                          effects = list(c(s.index, list(Intercept = 1)),
                                         list(dist = meuse$dist)),
                          tag = "meuse.data")

#Create data structure for prediction
A.pred <- inla.spde.make.A(mesh = mesh, loc = coordinates(meuse.grid))
meuse.stack.pred <- inla.stack(data = list(zinc = NA),
                               A = list(A.pred, 1),
                               effects = list(c(s.index, list (Intercept = 1)),
                                              list(dist = meuse.grid$dist)),
                               tag = "meuse.pred")

#Join stack
join.stack <- inla.stack(meuse.stack, meuse.stack.pred)

#Fit model
form <- log(zinc) ~ -1 + Intercept + dist + f(spatial.field, model = spde)

m1 <- inla(form, data = inla.stack.data(join.stack, spde = meuse.spde),
           family = "gaussian",
           control.predictor = list(A = inla.stack.A(join.stack), compute = TRUE),
           control.compute = list(cpo = TRUE, dic = TRUE))

#Summary of results
summary(m1)

#Get predicted data on grid
index.pred <- inla.stack.index(join.stack, "meuse.pred")$data

meuse.grid$zinc.spde <- m1$summary.fitted.values[index.pred, "mean"]
meuse.grid$zinc.spde.sd <- m1$summary.fitted.values[index.pred, "sd"]


#Compute statistics in terms or range and variance
spde.est <- inla.spde2.result(inla = m1, name = "spatial.field",
                              spde = meuse.spde, do.transf = TRUE)

#Kappa
#inla.zmarginal(spde.est$marginals.kappa[[1]])
#Variance
inla.zmarginal(spde.est$marginals.variance.nominal[[1]])

#Range
inla.zmarginal(spde.est$marginals.range.nominal[[1]])


