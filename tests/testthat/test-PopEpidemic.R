library(testthat)
library(AlphaSimR)
test_that("newPopEpidemic basic construction and defaults work", {

  set.seed(123)

  founderPop <- AlphaSimR::quickHaplo(nInd = 5, nChr = 1, segSites = 100)
  SP <- SimParamEpidemic$new(founderPop, model = "SIR")
  SP$addTraitA(nQtlPerChr = 10, mean = c(0,0,0), var = rep(1,3),
               corA = diag(rep(1,3)), name = c('sus', 'inf', 'tol'))
  
  ep <- newPopEpidemic(founderPop, simParam = SP)

  # Class inheritance
  expect_true(inherits(ep, "Pop"))

  # Dynamics data.table structure
  expect_true(data.table::is.data.table(ep@dynamics))
  # Required columns present (order not enforced)
  expect_true(all(c("indCases", "group", "status", "group_inf",
                    "generation", "infected_by") %in% names(ep@dynamics)))

  # Timings columns present and NA_real_
  if (length(SP$timings) > 0) {
    expect_true(all(SP$timings %in% names(ep@dynamics)))
    for (tcol in SP$timings) {
      expect_true(all(is.na(ep@dynamics[[tcol]])))
    }
  }

  # Defaults for grouping/indCases
  expect_all_equal(ep@dynamics$group, 1L)
  expect_equal(sum(ep@dynamics$indCases == 1), 1L)

  # Status factor levels and initial value
  expect_identical(levels(ep@dynamics$status), unique(strsplit("SIR", "")[[1]]))

  # group_inf default
  expect_true(all(ep@dynamics$group_inf == 0))

  # generation and infected_by values
  indCases <- which(ep@dynamics$indCases == 1)
  control_idx <- which(ep@dynamics$indCases == 0)
  expect_true(all(ep@dynamics$generation[indCases] == 1L))
  expect_true(all(ep@dynamics$infected_by[indCases] == 0L))
  expect_true(all(is.na(ep@dynamics$generation[control_idx])))
  expect_true(all(is.na(ep@dynamics$infected_by[control_idx])))
})

test_that("newPopEpidemic honors provided group_list and indCases", {

  founderPop <- AlphaSimR::quickHaplo(nInd = 5, nChr = 1, segSites = 50)
  SP <- SimParamEpidemic$new(founderPop, model = "SIR")
  SP$addTraitA(nQtlPerChr = 5, mean = c(0,0,0), var = c(1,1,1),
               cor = diag(c(1,1,1)), name = c('sus', 'inf', 'tol'))

  group_list <- c(1, 1, 2, 2, 2)
  indCases <- c(1, 0, 0, 1, 0)

  ep <- newPopEpidemic(founderPop, group_list = group_list,
                       indCases = indCases, simParam = SP)

  expect_identical(as.integer(ep@dynamics$group), as.integer(group_list))
  expect_identical(as.integer(ep@dynamics$indCases), as.integer(indCases))

  # One index per group
  expect_equal(sum(ep@dynamics$indCases[ep@dynamics$group == 1] == 1), 1L)
  expect_equal(sum(ep@dynamics$indCases[ep@dynamics$group == 2] == 1), 1L)

})

test_that("newPopEpidemic error handling works", {

  founderPop <- AlphaSimR::quickHaplo(nInd = 5, nChr = 1, segSites = 50)
  SPbase <- SimParam$new(founderPop)
  SP <- SimParamEpidemic$new(founderPop, model = "SIR")
  SP$addTraitA(nQtlPerChr = 5)

  # simParam must be SimParamEpidemic; AlphaSimR SimParam is not enough
  pop <- AlphaSimR::newPop(founderPop, simParam = SPbase)
  expect_error(asPopEpidemic(pop, simParam = SPbase),"SimParamEpidemic")
  

  # group_list length mismatch
  expect_error(newPopEpidemic(founderPop, group_list = rep(1L, 4), simParam = SP),
               "group_list must have equal elements")

  # indCases length mismatch
  expect_error(newPopEpidemic(founderPop,
                          group_list = rep(1L, 5),
                          indCases = c(1, 0, 0, 0),
                          simParam = SP),
               "indCases must have the same length")

  # indCases invalid values
  expect_error(newPopEpidemic(founderPop,
                          group_list = rep(1L, 5),
                          indCases = c(1, 0, 0, 2, 0),
                          simParam = SP),
               "indCases must contain only 0 and 1 values")

  # must have at least one index per group
  expect_error(newPopEpidemic(founderPop,
                          group_list = c(1, 1, 2, 2, 2),
                          indCases = c(1, 0, 0, 0, 0),
                          simParam = SP),
               "Each group must have at least one case")
})

test_that("show(PopEpidemic) prints expected summary", {
  skip_if_not_installed("AlphaSimR")

  founderPop <- AlphaSimR::quickHaplo(nInd = 3, nChr = 1, segSites = 50)
  SP <- SimParamEpidemic$new(founderPop, model = "SIR")
  ones <- rep(1.,3)
  SP$addTraitA(nQtlPerChr = 5, mean = ones-1, var = ones, corA = diag(ones),
               name= c('sus', 'inf', 'tol'))
  ep <- newPopEpidemic(founderPop, simParam = SP)
  expect_output(show(ep), "Model: SIR", fixed = TRUE)
})
