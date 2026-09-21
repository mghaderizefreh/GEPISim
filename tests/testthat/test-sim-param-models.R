test_that ("all supported models configure compartments, traits, and timings", {
  founders <- AlphaSimR::quickHaplo(nInd = 4, nChr = 1, segSites = 10)

  expected <- list(
    SIR = list(
      compartments = list("S", "I", "R"),
      epi_traits = c(s = "sus", i = "inf", t = "tol"),
      timings = c("Tinf", "Tdeath")
    ),
    SIDR = list(
      compartments = list("S", "I", "D", "R"),
      epi_traits = c(s = "sus", i = "inf", d = "det", t = "tol"),
      timings = c("Tinf", "Tsign", "Tdeath")
    ),
    SEIR = list(
      compartments = list("S", "E", "I", "R"),
      epi_traits = c(s = "sus", i = "inf", l = "lat", t = "tol"),
      timings = c("Tinf", "Tsign", "Tdeath")
    ),
    SEIDR = list(
      compartments = list("S", "E", "I", "D", "R"),
      epi_traits = c(s = "sus", i = "inf", l = "lat", d = "det", t = "tol"),
      timings = c("Tinf", "Tinc", "Tsign", "Tdeath")
    )
  )

  for (model in names(expected)) {
    sp <- SimParamEpidemic$new(founders, model = model)

    expect_identical(sp$model, model)
    expect_identical(sp$compartments, expected[[model]]$compartments)
    expect_identical(sp$epi_traits, expected[[model]]$epi_traits)
    expect_identical(sp$timings, expected[[model]]$timings)
  }
})

test_that("version is a read-only active field", {
  founders <- AlphaSimR::quickHaplo(nInd = 4, nChr = 1, segSites = 10)
  sp <- SimParamEpidemic$new(founders)

  expect_named(sp$version, c("AlphaSimR", "GEPISim"))
  expect_error(
    sp$version <- list(),
    "`\\$version` is read only"
  )
})

