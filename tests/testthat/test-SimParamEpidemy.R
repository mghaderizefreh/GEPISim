library(testthat)
library(AlphaSimR)

# tests from AlphaSimR applied to epiAlphaSimR ----
test_that("SimParamEpidemy nThreads validates values and NULL resets to default", {
  founder <- quickHaplo(nInd = 2, nChr = 2, segSites = 4)
  SP <- SimParamEpidemy$new(founder)
  pop <- newEpidemy(founder, simParam = SP)

  SP$nThreads <- 1L
  expect_equal(SP$nThreads, 1L)

  SP$nThreads <- NULL
  expect_equal(SP$nThreads, getNumThreads())
  expect_silent(pullSegSiteGeno(pop, simParam = SP))

  expect_error(
    SP$nThreads <- 0,
    regexp = "single positive integer or NULL to reset"
  )
  expect_error(
    SP$nThreads <- 0L,
    regexp = "single positive integer or NULL to reset"
  )
  expect_error(
    SP$nThreads <- 1.5,
    regexp = "single positive integer or NULL to reset"
  )
  expect_error(
    SP$nThreads <- NA_integer_,
    regexp = "single positive integer or NULL to reset"
  )
})

#### new tests for epiAlphaSimR ----
make_founders <- function() {
  # Small/fast founder genomes for tests
  AlphaSimR::quickHaplo(nInd = 4, nChr = 1, segSites = 10)
}

test_that("constructor creates an object of the correct class and inheritance", {
  founders <- make_founders()
  sp <- SimParamEpidemy$new(founders)

  expect_true(inherits(sp, "SimParamEpidemy"))
  expect_true(inherits(sp, "SimParam"))  # inherits from AlphaSimR::SimParam (R6)
})

test_that("constructor sets default arguments correctly", {
  founders <- make_founders()
  sp <- SimParamEpidemy$new(founders)

  expect_identical(sp$model, "SIR")
  expect_identical(sp$removal_period, 10)
  expect_identical(sp$r_beta, 0.5)
  expect_identical(sp$RP_shape, 1)
  expect_identical(sp$RP_scale, 10)  # removal_period / RP_shape

  # Model-dependent defaults for SIR
  expect_identical(sp$epi_traits, c(s = "sus", i = "inf", t = "tol"))
  expect_identical(sp$timings, c("Tinf", "Tdeath"))
})

test_that("constructor sets custom arguments correctly", {
  founders <- make_founders()
  model <- "SIR"
  removal_period <- 8
  r_beta <- 0.75
  RP_shape <- 2

  sp <- SimParamEpidemy$new(
    founders,
    model = model,
    removal_period = removal_period,
    r_beta = r_beta,
    RP_shape = RP_shape
  )

  expect_identical(sp$model, model)
  expect_identical(sp$removal_period, removal_period)
  expect_identical(sp$r_beta, r_beta)
  expect_identical(sp$RP_shape, RP_shape)
  expect_identical(sp$RP_scale, removal_period / RP_shape)

  # Model-dependent defaults for SIR still apply
  expect_identical(sp$epi_traits, c(s = "sus", i = "inf", t = "tol"))
  expect_identical(sp$timings, c("Tinf", "Tdeath"))
})

test_that("constructor rejects invalid model values", {
  founders <- make_founders()

  expect_error(
    SimParamEpidemy$new(founders, model = "ABC"),
    regexp = "provided model is not valid"
  )
})

test_that("constructor accepts model value case-insensitively (as documented)", {
  founders <- make_founders()

  # As documented, model is case-insensitive.
  # This test will currently fail with the present implementation and should
  # pass once validation is adjusted to be case-insensitive.
  expect_silent({
    sp <- SimParamEpidemy$new(founders, model = "sir")
    # Ensure model-dependent components are set correctly even if lower-case was passed
    expect_identical(sp$epi_traits, c(s = "sus", i = "inf", t = "tol"))
    expect_identical(sp$timings, c("Tinf", "Tdeath"))
  })
})

test_that("active binding 'version' returns expected structure and values", {
  founders <- make_founders()
  sp <- SimParamEpidemy$new(founders)

  v <- sp$version
  expect_type(v, "list")
  expect_true(all(c("AlphaSimR", "EpiAlphaSimR") %in% names(v)))

  # EpiAlphaSimR version should match installed package version
  expect_identical(
    v$EpiAlphaSimR,
    utils::packageDescription("epiAlphaSimR")$Version
  )
})
