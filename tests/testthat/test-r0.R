test_that ("get_R0 returns the median individual R0 for SIR", {
  x <- make_test_epidemic(
    model = "SIR",
    n_ind = 3,
    ind_cases = c(1L, 0L, 0L)
  )

  x$sp$r_beta <- 0.1
  x$sp$removal_period <- 10

  x$epop@pheno[, "sus"] <- c(1, 2, 4)
  x$epop@pheno[, "inf"] <- c(2, 3, 5)
  x$epop@pheno[, "tol"] <- c(1, 0.5, 2)

  expected_r0s <- vapply(seq_len(x$epop@nInd), function(i) {
    x$sp$r_beta *
      mean(x$epop@pheno[-i, "sus"]) *
      x$epop@pheno[i, "inf"] *
      x$sp$removal_period *
      x$epop@pheno[i, "tol"]
  }, numeric(1))

  expect_identical(
    get_R0(x$epop, simParam = x$sp),
    median(expected_r0s)
  )
})

test_that("get_R0 returns the median individual R0 for SIDR", {
  x <- make_test_epidemic(
    model = "SIDR",
    n_ind = 3,
    ind_cases = c(1L, 0L, 0L)
  )
  
  x$sp$r_beta <- 0.1
  x$sp$removal_period <- 10
  x$sp$detection_period <- 4
  
  x$epop@pheno[, "sus"] <- c(1, 2, 4)
  x$epop@pheno[, "inf"] <- c(2, 3, 5)
  x$epop@pheno[, "tol"] <- c(1, 0.5, 2)
  x$epop@pheno[, "det"] <- c(0.5, 2, 1.5)
  
  # The SIDR implementation defines individual R0 as:
  #
  # beta * mean(susceptibility of everyone except i) * infectivity_i *
  #   (removal_period * tolerance_i + detection_period * detection_i)
  #
  # Individual 1:
  # 0.1 * mean(2, 4) * 2 * (10 * 1 + 4 * 0.5) = 7.2
  #
  # Individual 2:
  # 0.1 * mean(1, 4) * 3 * (10 * 0.5 + 4 * 2) = 9.75
  #
  # Individual 3:
  # 0.1 * mean(1, 2) * 5 * (10 * 2 + 4 * 1.5) = 19.5
  #
  # Therefore median(c(7.2, 9.75, 19.5)) = 9.75.
  expected_r0 <- 9.75
  
  expect_equal(
    get_R0(x$epop, simParam = x$sp),
    expected_r0
  )
})

