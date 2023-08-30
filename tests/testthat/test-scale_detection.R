library(SpatialExperiment)

test_that("max octave works", {
  spe = SpatialExperiment(
    assays = list("counts" = matrix(0, 2, 2)),
    colData = data.frame(x = c(1, 500), y = c(1, 500)),
    spatialCoordsNames = c("x", "y")
  )

  expect_equal(getOctaves(spe), 5)
  spatialCoords(spe)[2, ] = 10000
  expect_equal(getOctaves(spe), 9)
  
  # multiple sample ids
  spe = SpatialExperiment(
    assays = list("counts" = matrix(0, 2, 2)),
    colData = data.frame(x = c(1, 500), y = c(1, 500), sample_id = c('s1', 's2')),
    spatialCoordsNames = c("x", "y")
  )
  expect_error(getOctaves(spe))
})

test_that("max scale works", {
  expect_equal(getMaxOctaveScales(1), 27)
  expect_equal(getMaxOctaveScales(2), 51)
  expect_equal(getMaxOctaveScales(3), 99)
  expect_equal(getMaxOctaveScales(4), 195)
})

test_that("mask check works", {
  spe = SpatialExperiment(
    assays = list("counts" = matrix(0, 2, 2)),
    colData = data.frame(x = c(1, 5), y = c(1, 10)),
    spatialCoordsNames = c("x", "y")
  )

  mask = matrix(rnorm(50), 10, 5)
  expect_error(checkMask(mask, spe), "binary")
  expect_silent(checkMask(mask > 0, spe))
  expect_silent(checkMask(apply(mask > 0, 2, as.numeric), spe))
  mask[1, 1] = NA_real_
  expect_error(checkMask(mask > 0, spe), "missing")
  expect_error(checkMask(cbind(mask, 1), spe), "size")
  expect_error(checkMask(rbind(mask, 1), spe), "size")

})
