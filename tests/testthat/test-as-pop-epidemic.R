test_that ("default index-case selection creates one index case per group", {
  founders <- AlphaSimR::quickHaplo(nInd = 6, nChr = 1, segSites = 20)

  sp <- SimParamEpidemic$new(founders, model = "SIR")
  sp$addTraitA(
    nQtlPerChr = 2,
    mean = c(0, 0, 0),
    var = c(1, 1, 1),
    corA = diag(3),
    name = c("sus", "inf", "tol")
  )

  group_list <- c("herd-A", "herd-A", "herd-B", "herd-B", "herd-C", "herd-C")

  set.seed(100)
  epop <- newPopEpidemic(
    founders,
    group_list = group_list,
    simParam = sp
  )

  index_cases_per_group <- epop@dynamics[
    ,
    sum(indCases == 1L),
    by = group
  ]

  expect_identical(index_cases_per_group$group, unique(group_list))
  expect_true(all(index_cases_per_group$V1 == 1L))
  expect_equal(sum(epop@dynamics$indCases == 1L), 3L)
})

test_that("asPopEpidemic rejects populations missing required epidemic phenotypes", {
  founders <- AlphaSimR::quickHaplo(nInd = 4, nChr = 1, segSites = 20)

  sp <- SimParamEpidemic$new(founders, model = "SIR")

  # Deliberately omit the phenotype named "tol", which SIR requires.
  sp$addTraitA(
    nQtlPerChr = 2,
    mean = c(0, 0, 0),
    var = c(1, 1, 1),
    corA = diag(3),
    name = c("sus", "inf", "not_tol")
  )

  pop <- AlphaSimR::newPop(founders, simParam = sp)

  expect_error(
    asPopEpidemic(pop, simParam = sp),
    '"tol" is not defined as a phenotype'
  )
})

test_that("initial dynamics use model-specific status levels", {
  models <- c("SIR", "SIDR", "SEIR", "SEIDR")

  for (model in models) {
    x <- make_test_epidemic(model = model)

    expect_identical(
      levels(x$epop@dynamics$status),
      unlist(x$sp$compartments)
    )
    expect_true(all(as.character(x$epop@dynamics$status) == "S"))
    expect_true(all(x$sp$timings %in% names(x$epop@dynamics)))
  }
})

