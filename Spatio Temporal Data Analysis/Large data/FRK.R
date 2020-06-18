library('sp')
library(FRK)
library(dggrids)
library(dplyr)

# devtools::install_github("hadley/devtools")
# devtools::install_github("andrewzm/dggrids")

tol = 0.01
f = TrueTemp ~ 1
BAUs <- all.sat.temps                    # assign BAUs
BAUs$Missing <- is.na(BAUs$MaskTemp)      # mark which BAUs contain missing data
BAUs$MaskTemp <- NULL                     # remove data from BAU
BAUs$TrueTemp <- NULL
BAUs$fs <- 1                          # set fs variation to unity
coordinates(BAUs)  <- ~Lon+Lat        # convert to SpatialPointsDataFrame
gridded(BAUs) <- TRUE  


train.FRK <- subset(all.sat.temps,!is.na(MaskTemp))        # no missing data in data frame
coordinates(train.FRK)  <- ~Lon+Lat

basis <- auto_basis(plane(),          # we are on the plane
                    data = train.FRK,       # data around which to make basis
                    regular = 0,      # irregular basis
                    nres = 3,         # 3 resolutions
                    scale_aperture = 1)   # aperture scaling of basis functions 
S <- FRK(f = f,                       # formula for SRE model
         data = train.FRK,                  # data
         basis = basis,               # Basis
         BAUs = BAUs,                 # BAUs
         tol = tol) 
S

## Predict
BAUs_pred <- predict(S)           # predict over all BAUs
BAUs_pred_df <- data.frame(BAUs_pred) # convert to data frame
tmp = BAUs_pred_df %>% filter(Missing == T)

sat.test$frk_fitted = tmp$mu
sat.test$frk_fitted_sd = tmp$sd

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

evaluate(sat.test$TrueTemp, sat.test$frk_fitted, sat.test$frk_fitted_sd)

save.image('final_workspace.RData')
