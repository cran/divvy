test_that('rangeSize returns appropriate values', {
  xy <- cbind('x' = c(1, 1, 1),
              'y' = c(0, 5, 10))
  rng <- rangeSize(xy)
  expect_equal(rng[, 1:4], c(3, 1, 5, 10), ignore_attr = TRUE)
  expect_type(rng[, 5:7], 'double')
  xy <- data.frame(xy)
  rng2 <- rangeSize(xy)
  expect_identical(rng, rng2)
})

test_that('return NA distances for single pt', {
  pt <- cbind(1, 1)
  rngPt <- rangeSize(pt)
  expect_equal(sum(is.na(rngPt)), 4)
})

test_that('sdSumry returns values within basic tolerance limits', {
  smry <- sdSumry(bivalves, c('paleolng','paleolat'), 'genus', quotaQ = 0.8, quotaN = 100)
  expect_s3_class(smry, 'data.frame')
  expect_gt(smry[, 'SQSdiv'], 200)
  expect_gt(smry[, 'CRdiv'],  300)
  expect_gt(smry[, 'SQSupr95'], smry[, 'SQSlow95'])
  expect_gt(smry[, 'CRupr95'],  smry[, 'CRlow95'])
  expect_equal(smry[, 'nOcc'], 5096)
  expect_equal(smry[, 'meanPairDist'], 7277.259, tolerance = 0.01)
})

# project occurrences and retrieve cell centroids in new coordinate system
rWorld <- terra::rast()
prj <- 'EPSG:8857'
rPrj <- terra::project(rWorld, prj, res = 200000)
xyCartes <- c('paleolng','paleolat')
xyCell   <- c('centroidX','centroidY')
llOccs <- terra::vect(bivalves, geom = xyCartes, crs = 'epsg:4326')
prjOccs <- terra::project(llOccs, prj)
cellIds <- terra::cells(rPrj, prjOccs)[,'cell']
bivAlt <- bivalves
bivAlt[, xyCell] <- terra::xyFromCell(rPrj, cellIds)

test_that('rangeSize and sdSumry accept projected coords', {
  expect_no_condition( rangeSize(bivAlt[, xyCell], crs = prj) )
  rngPrj <- rangeSize(bivAlt[, xyCell], crs = prj)
  expect_type(rngPrj, 'double')
  expect_equal(sum(is.na(rngPrj)), 0)

  expect_no_condition( sdSumry(bivAlt,   xyCell,  'genus', crs = prj) )
  smryLL  <- sdSumry(bivalves, xyCartes, 'genus')
  smryPrj <- sdSumry(bivAlt,   xyCell,   'genus', crs = prj)
  expect_equal(smryLL[,'latRange'], smryPrj[,'latRange'], tolerance = 0.01)
})
