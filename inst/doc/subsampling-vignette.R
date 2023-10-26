## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = FALSE,
  comment = "#>"
)

## ----setup, message=FALSE-----------------------------------------------------
library(divvy) 
library(units)
library(sf)
library(terra)
library(vegan)
library(iNEXT)

## ----attached data------------------------------------------------------------
data(bivalves)
head(bivalves)
nrow(bivalves)
length(unique(bivalves$genus))

## ----rasterise----------------------------------------------------------------
# initialise Equal Earth projected coordinates
rWorld <- rast()
prj <- 'EPSG:8857'
rPrj <- project(rWorld, prj, res = 200000) # 200,000m is approximately 2 degrees
values(rPrj) <- 1:ncell(rPrj)

# coordinate column names for the current and target coordinate reference system
xyCartes <- c('paleolng','paleolat')
xyCell   <- c('cellX','cellY')

# retrieve coordinates of raster cell centroids
llOccs <- vect(bivalves, geom = xyCartes, crs = 'epsg:4326')
prjOccs <- project(llOccs, prj)
bivalves$cell <- cells(rPrj, prjOccs)[,'cell']
bivalves[, xyCell] <- xyFromCell(rPrj, bivalves$cell)

## ----map data, fig.width=6, fig.align='center', message=FALSE-----------------
occUniq <- uniqify(bivalves, xyCell)
ptsUniq <- st_as_sf(occUniq, coords = xyCell, crs = prj)

library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2) # for plot visualisation
world <- ne_countries(scale = "medium", returnclass = "sf")
worldP <- ggplot(data = world) +
  theme_bw() +
  geom_sf() +
  geom_sf(data = ptsUniq, shape = 17, color = 'blue')

plot(worldP)

## ----circ subsample-----------------------------------------------------------
set.seed(1)
circLocs <- cookies(dat = bivalves,  
                    xy = xyCell,
                    iter = 500, 
                    nSite = 12, 
                    r = 1500, # radial distance in km
                    weight = TRUE, # probabilistically aggregate subsampling sites
                    crs = prj, # Equal Earth projection
                    output = 'locs')
length(circLocs)
circLocs[[1]]

## ----map subsample, fig.width=6, fig.align='center'---------------------------
# over-plot the subsample locations
smplPts <- st_as_sf(circLocs[[1]], coords = xyCell, crs = prj)
worldP +
  geom_sf(data = smplPts, shape = 17, color = 'red')

## ----count pool size, fig.width=6, fig.align='center'-------------------------
cntr <- smplPts[1,]
# distances inferred to be meters based on lat-long coord system
r <- 1500
buf <- st_buffer(cntr, dist = r*1000)

# over-plot subsample boundary region and included sites
worldP +
  geom_sf(data = smplPts, shape = 17, color = 'red') +
  geom_sf(data = buf, fill = NA, linewidth = 1, color = 'red')

## ----tally pool size----------------------------------------------------------
inBuf <- st_intersects(ptsUniq, buf, sparse = FALSE) 
sum(inBuf)

## ----circ subsample variation-------------------------------------------------
set.seed(7)
# same parameter values as above except for 'output'
circOccs <- cookies(dat = bivalves, 
                    xy = xyCell, 
                    iter = 500, 
                    nSite = 12, 
                    r = 1500, 
                    weight = TRUE, 
                    crs = prj,
                    output = 'full')
head( circOccs[[8]] )

## ----unique occs--------------------------------------------------------------
bivUniq <- uniqify(bivalves, taxVar = 'genus', xyCell)
nrow(bivUniq)

## ----MST subsample------------------------------------------------------------
set.seed(2)
nnLocs <- clustr(dat = bivalves, 
                 xy = xyCell, 
                 iter = 500, 
                 distMax = 3000, # diameter = 2x the circular radius set above
                 nSite = 12,
                 crs = prj
                 )
nnLocs[[1]]

## ----MST all sites------------------------------------------------------------
set.seed(2)
nnAllSites <- clustr(dat = bivalves, 
                     xy = xyCell, 
                     iter = 500,
                     distMax = 3000, # diameter = 2x the circular radius set above
                     nMin = 3,
                     crs = prj
                     )
nrow( nnAllSites[[1]] )

## ----lat subsample------------------------------------------------------------
bandLocs <- bandit(dat = bivalves,
                   xy = xyCell,
                   iter = 100, nSite = 12, 
                   bin = 20, # interval width in degrees
                   crs = prj
                   # ,absLat = TRUE # optional
                   )
nrow(bandLocs[[1]]) # number of sites in one subsampled band
length(bandLocs) # number of subsamples, tallied across all bands
unique(names(bandLocs)) # latitudinal degrees of subsampled intervals

## ----meta data----------------------------------------------------------------
unsamp <- sdSumry(dat = bivalves, 
                  taxVar = 'genus',
                  collections = 'collection_no',
                  xy = xyCell,
                  quotaQ = 0.8, quotaN = 100,
                  omitDom = TRUE,
                  crs = prj)
unsamp

samp1 <- sdSumry(dat = circOccs[[1]], 
                 taxVar = 'genus',
                 collections = 'collection_no',
                 xy = xyCell,
                 quotaQ = 0.8, quotaN = 100, 
                 omitDom = TRUE,
                 crs = prj)
samp1

## ----summary over subsamples--------------------------------------------------
# warning - it's slow to rarefy diversity for hundreds of subsamples!
# this code chunk skips it for quick demonstration purposes
sampsMeta <- sdSumry(dat = circOccs, 
                     taxVar = 'genus',
                     collections = 'collection_no',
                     xy = xyCell,
                     crs = prj
                     )
quantile(sampsMeta$nTax, c(0.25, 0.5, 0.75))

