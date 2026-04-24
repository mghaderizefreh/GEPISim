library(testthat)
library(AlphaSimR)
test_that("newEpidemy basic construction and defaults work", {

  set.seed(123)

  founderPop <- AlphaSimR::quickHaplo(nInd = 5, nChr = 1, segSites = 100)
  SP <- SimParamEpidemy$new(founderPop, model = "SIR")
  SP$addTraitA(nQtlPerChr = 10)

  ep <- newEpidemy(founderPop, simParam = SP)

  # Class and basic slots
  expect_s4_class(ep, "EPop")
  expect_true(inherits(ep, "Pop"))
  expect_equal(ep@nInd, founderPop@nInd)
  expect_equal(ep@nTraits, SP$nTraits)

  # Dynamics data.table structure
  expect_true(data.table::is.data.table(ep@dynamics))
  # Required columns present (order not enforced)
  expect_true(all(c("donor", "group", "status", "group_inf",
                    "generation", "infected_by") %in% names(ep@dynamics)))

  # Timings columns present and NA_real_
  if (length(SP$timings) > 0) {
    expect_true(all(SP$timings %in% names(ep@dynamics)))
    for (tcol in SP$timings) {
      expect_true(all(is.na(ep@dynamics[[tcol]])))
      # NAs should be of real type
      expect_true(is.double(ep@dynamics[[tcol]]) || all(is.na(ep@dynamics[[tcol]])))
    }
  }

  # Defaults for grouping/donors
  expect_all_equal(ep@dynamics$group, 1L)
  expect_equal(sum(ep@dynamics$donor == 1), 1L)

  # Status factor levels and initial value
  expect_identical(levels(ep@dynamics$status), unique(strsplit("SIR", "")[[1]]))

  # TODO: This actually should not be true as the donor is not susceptible but
  #       may depend on the model if the model is infectious only or also
  #       diseased as well. Need to confirm with Jamie or see his code
  #  expect_true(all(as.character(ep@dynamics$status) == "S"))

  # group_inf default
  expect_true(all(ep@dynamics$group_inf == 0))

  # generation and infected_by values
  donors_idx <- which(ep@dynamics$donor == 1)
  nondonors_idx <- which(ep@dynamics$donor == 0)
  expect_true(all(ep@dynamics$generation[donors_idx] == 1L))
  expect_true(all(ep@dynamics$infected_by[donors_idx] == 0L))
  expect_true(all(is.na(ep@dynamics$generation[nondonors_idx])))
  expect_true(all(is.na(ep@dynamics$infected_by[nondonors_idx])))
})

test_that("newEpidemy honors provided group_list and donor_list", {

  founderPop <- AlphaSimR::quickHaplo(nInd = 5, nChr = 1, segSites = 50)
  SP <- SimParamEpidemy$new(founderPop, model = "SIR")
  SP$addTraitA(nQtlPerChr = 5)

  group_list <- c(1, 1, 2, 2, 2)
  donor_list <- c(1, 0, 0, 1, 0)

  ep <- newEpidemy(founderPop, group_list = group_list,
                   donor_list = donor_list, simParam = SP)

  expect_identical(as.integer(ep@dynamics$group), as.integer(group_list))
  expect_identical(as.integer(ep@dynamics$donor), as.integer(donor_list))

  # One donor per group
  expect_equal(sum(ep@dynamics$donor[ep@dynamics$group == 1] == 1), 1L)
  expect_equal(sum(ep@dynamics$donor[ep@dynamics$group == 2] == 1), 1L)

  # generation/infected_by at donors vs non-donors
  donors_idx <- which(ep@dynamics$donor == 1)
  nondonors_idx <- which(ep@dynamics$donor == 0)
  expect_true(all(ep@dynamics$generation[donors_idx] == 1L))
  expect_true(all(ep@dynamics$infected_by[donors_idx] == 0L))
  expect_true(all(is.na(ep@dynamics$generation[nondonors_idx])))
  expect_true(all(is.na(ep@dynamics$infected_by[nondonors_idx])))
})

test_that("newEpidemy error handling works", {

  founderPop <- AlphaSimR::quickHaplo(nInd = 5, nChr = 1, segSites = 50)
  SPbase <- SimParam$new(founderPop)
  SP <- SimParamEpidemy$new(founderPop, model = "SIR")
  SP$addTraitA(nQtlPerChr = 5)

  # simParam must be SimParamEpidemy; AlphaSimR SimParam is not enough
  expect_error(newEpidemy(founderPop, simParam = SPbase),
               "SimParamEpidemy")

  # group_list length mismatch
  expect_error(newEpidemy(founderPop, group_list = rep(1L, 4), simParam = SP),
               "group_list must have equal elements")

  # donor_list length mismatch
  expect_error(newEpidemy(founderPop,
                          group_list = rep(1L, 5),
                          donor_list = c(1, 0, 0, 0),
                          simParam = SP),
               "donor_list must have the same length")

  # donor_list invalid values
  expect_error(newEpidemy(founderPop,
                          group_list = rep(1L, 5),
                          donor_list = c(1, 0, 0, 2, 0),
                          simParam = SP),
               "donor_list must contain only 0 and 1 values")

  # must have at least one donor per group
  expect_error(newEpidemy(founderPop,
                          group_list = c(1, 1, 2, 2, 2),
                          donor_list = c(1, 0, 0, 0, 0),
                          simParam = SP),
               "Each group must have at least one donor")
})

test_that("show(EPop) prints expected summary", {
  skip_if_not_installed("AlphaSimR")

  founderPop <- AlphaSimR::quickHaplo(nInd = 3, nChr = 1, segSites = 50)
  SP <- SimParamEpidemy$new(founderPop, model = "SIR")
  SP$addTraitA(nQtlPerChr = 5)
  ep <- newEpidemy(founderPop, simParam = SP)

  expect_output(show(ep), "EPop")
  expect_output(show(ep), "Ploidy:")
  expect_output(show(ep), "Individuals:")
  expect_output(show(ep), "Chromosomes:")
  expect_output(show(ep), "Loci:")
  expect_output(show(ep), "Traits:")
  expect_output(show(ep), "Model: SIR")
})