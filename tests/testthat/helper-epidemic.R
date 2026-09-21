make_test_epidemic <- function(model = "SIR",
                               n_ind = 6L,
                               group_list = rep("group-1", n_ind),
                               ind_cases = c(1L, rep(0L, n_ind - 1L)),
                               r_beta = 0) {
  stopifnot(length(group_list) == n_ind, length(ind_cases) == n_ind)

  founders <- AlphaSimR::quickHaplo(
    nInd = n_ind,
    nChr = 1,
    segSites = 20
  )

  sp <- SimParamEpidemic$new(
    founders,
    model = model,
    r_beta = r_beta,
    removal_period = 2,
    detection_period = 2,
    latent_period = 2
  )

  trait_names <- unname(sp$epi_traits)
  n_traits <- length(trait_names)

  sp$addTraitA(
    nQtlPerChr = 2,
    mean = rep(0, n_traits),
    var = rep(1, n_traits),
    corA = diag(n_traits),
    name = trait_names
  )

  epop <- newPopEpidemic(
    founders,
    group_list = group_list,
    indCases = ind_cases,
    simParam = sp
  )

  # The epidemic models require strictly positive values because these values
  # scale gamma-distributed waiting times and infection rates.
  epop@pheno[, trait_names] <- 1

  list(epop = epop, sp = sp)
}

