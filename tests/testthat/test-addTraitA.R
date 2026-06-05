library(testthat)
library(AlphaSimR)

## tests from AlphaSimR package ----
#Population with 2 individuals, 1 chromosome and 1 QTL
#Population is fully inbred and p=q=0.5
founderPop = newMapPop(list(c(0)),
                       list(matrix(c(1,1,0,0),
                                   nrow=4,ncol=1)))

test_that("addTraitA",{
  SP = SimParamEpidemic$new(founderPop=founderPop, model = "SIR")
  SP$nThreads = 1L
  SP$addTraitA(nQtlPerChr=1,mean=0,var=1)
  pop = newPop(founderPop,simParam=SP)
  expect_equal(abs(SP$traits[[1]]@addEff),1,tolerance=1e-6)
  expect_equal(abs(SP$traits[[1]]@intercept),0,tolerance=1e-6)
  expect_equal(SP$varA,rep(1,nchar(SP$model)),tolerance=1e-6) # adapted
  expect_equal(SP$varG,rep(1,nchar(SP$model)),tolerance=1e-6) # adapted
  ans = genParam(pop,simParam=SP)
  expect_equal(abs(unname(c(ans$varA))),
               rep(1,nchar(SP$model)^2),tolerance=1e-6) # adapted
  expect_equal(abs(unname(c(ans$varD))),
               rep(0,nchar(SP$model)^2),tolerance=1e-6) # adapted
  expect_equal(abs(unname(c(ans$varG))),
               rep(1,nchar(SP$model)^2),tolerance=1e-6) # adapted
  expect_equal(unname(ans$genicVarA),rep(0.5,nchar(SP$model)),tolerance=1e-6)#ad
  expect_equal(unname(ans$genicVarD),rep(0,nchar(SP$model)),tolerance=1e-6) #ad
  expect_equal(unname(ans$genicVarG),rep(0.5,nchar(SP$model)),tolerance=1e-6)#ad
})

## New tests ----

# Helper to create a fresh SimParamEpidemic with a small founder population
.new_SP <- function(model = "SIR") {
  founderPop <- AlphaSimR::quickHaplo(nInd = 5, nChr = 1, segSites = 10)
  SP <- SimParamEpidemic$new(founderPop, model = model)
  SP$nThreads <- 1L
  SP
}

# Helper to extract the allowable trait "names" (values) from the model
.epi_trait_values <- function(SP) {
  unname(SP$epi_traits)
}

test_that("addTraitA: defaults to epi_traits when name=NULL and returns invisibly", {
  SP <- .new_SP()
  expect_invisible(SP$addTraitA(nQtlPerChr = 1))
})

test_that("addTraitA: accepts explicit names matching epi_traits (full set)", {
  SP <- .new_SP()
  trait_vals <- .epi_trait_values(SP)
  expect_silent(
    SP$addTraitA(
      nQtlPerChr = 1,
      name = trait_vals,
      mean = 0,
      var = 1,
      corA = diag(length(trait_vals))
    )
  )
})

test_that("addTraitA: rejects names that do not match epi_traits values", {
  SP <- .new_SP()
  # 'foo', 'bar' and 'baz' are not in the allowed epi_traits values
  expect_error(
    SP$addTraitA(nQtlPerChr = 1, name = c("foo", "bar", "baz")),
    "Provided names do not match the names from the model"
  )
})

test_that("addTraitA: rejects short-code names like 's','i','t' (must use epi_traits values)", {
  SP <- .new_SP()
  # For SIR, epi_traits values are  c("sus","inf","tol"), not "s","i","t"
  expect_error(
    SP$addTraitA(nQtlPerChr = 1, name = c("s", "i", "t")),
    "Provided names do not match the names from the model"
  )
})

test_that("addTraitA: accepts names that match epi_traits values even if they are not in the typical order", {
  SP <- .new_SP()
  # For SIR, epi_traits values are  c("sus","inf","tol"), but we provide them in a different order
  expect_silent(
    SP$addTraitA(nQtlPerChr = 1, name = c("tol", "sus", "inf"), mean = 0, var = 1)
  )
})

test_that("addTraitA: warns and forces means to zero when non-zero means are passed", {
  SP <- .new_SP()
  # Use a single trait for simplicity
  one_name <- .epi_trait_values(SP)[1]
  expect_warning(
    SP$addTraitA(nQtlPerChr = 1, name = one_name, mean = 1, var = 1),
    "Non zero mean values were passed. They will be forced to zero"
  )
})

test_that("addTraitA: var must be scalar or have length equal to length(name)", {
  SP <- .new_SP()
  trait_vals <- .epi_trait_values(SP)
  # Length mismatch: 2 variances for 3 names (for SIR)
  expect_error(
    SP$addTraitA(nQtlPerChr = 1, name = trait_vals, var = c(1, 1))
  )
  # Accepts scalar var for multiple names
  expect_silent(
    SP$addTraitA(nQtlPerChr = 1, name = trait_vals, var = 2)
  )
  # Accepts var vector with correct length
  SP2 <- .new_SP()
  expect_silent(
    SP2$addTraitA(nQtlPerChr = 1, name = trait_vals, var = rep(3, length(trait_vals)))
  )
})

test_that("addTraitA: corA must be NULL or a square matrix of size length(name)", {
  SP <- .new_SP()
  trait_vals <- .epi_trait_values(SP)
  # Wrong-sized corA (2x2 for 3 traits)
  expect_error(
    SP$addTraitA(nQtlPerChr = 1, name = trait_vals, corA = diag(2)),
    "corA should be NULL or square matrix of size length"
  )
  # Correct-sized corA
  SP2 <- .new_SP()
  expect_silent(
    SP2$addTraitA(nQtlPerChr = 1, name = trait_vals, corA = diag(length(trait_vals)))
  )
  # corA = NULL should be accepted (internally becomes identity)
  SP3 <- .new_SP()
  expect_silent(
    SP3$addTraitA(nQtlPerChr = 1, name = trait_vals, corA = NULL)
  )
})

test_that("addTraitA: scalar mean and var are recycled to length(name) before validation", {
  SP <- .new_SP()
  trait_vals <- .epi_trait_values(SP)
  # mean = 0 is silently recycled to the appropriate length
  expect_silent(
    SP$addTraitA(nQtlPerChr = 1, name = trait_vals, mean = 0, var = 1)
  )
})