test_that(
  "runEpidemic dispatches correctly and preserves no-transmission outbreaks", {
    expected_timings <- list(
      SIR = c("Tinf", "Tdeath"),
      SIDR = c("Tinf", "Tsign", "Tdeath"),
      SEIR = c("Tinf", "Tsign", "Tdeath"),
      SEIDR = c("Tinf", "Tinc", "Tsign", "Tdeath")
    )
    
    expected_final_status <- list(
      SIR = "R",
      SIDR = "R",
      SEIR = "R",
      SEIDR = "R"
    )
    
    group_list <- rep(c("group-A", "group-B"), each = 3)
    ind_cases <- c(1L, 0L, 0L, 1L, 0L, 0L)
    
    for (model in names(expected_timings)) {
      x <- make_test_epidemic(
        model = model,
        group_list = group_list,
        ind_cases = ind_cases,
        r_beta = 0
      )
      
      result <- suppressMessages(runEpidemic(x$epop, simParam = x$sp))
      dynamics <- result@dynamics
      
      index_rows <- which(dynamics$indCases == 1L)
      susceptible_rows <- which(dynamics$indCases == 0L)
      
      expect_s4_class(result, "PopEpidemic")
      expect_identical(result@iid, x$epop@iid)
      expect_identical(result@dynamics$group, x$epop@dynamics$group)
      
      # With beta = 0, only supplied index cases can ever be infected.
      expect_true(all(as.character(dynamics$status[index_rows]) ==
                        expected_final_status[[model]]))
      expect_true(all(as.character(dynamics$status[susceptible_rows]) == "S"))
      
      expect_true(all(dynamics$generation[index_rows] == 1L))
      expect_true(all(dynamics$infected_by[index_rows] == 0L))
      expect_true(all(is.na(dynamics$generation[susceptible_rows])))
      expect_true(all(is.na(dynamics$infected_by[susceptible_rows])))
      
      # Every seeded case obtains all model-specific event times.
      for (timing in expected_timings[[model]]) {
        expect_true(all(is.finite(dynamics[[timing]][index_rows])))
        expect_true(all(is.na(dynamics[[timing]][susceptible_rows])))
      }
      
      # Infection starts at time zero for supplied index cases.
      expect_identical(dynamics$Tinf[index_rows], rep(0, length(index_rows)))
    }
  })

test_that ("event times follow each model's disease progression order", {
  group_list <- rep(c("group-A", "group-B"), each = 3)
  ind_cases <- c(1L, 0L, 0L, 1L, 0L, 0L)
  
  # SIR: infection -> removal
  x <- make_test_epidemic("SIR", group_list = group_list,
                          ind_cases = ind_cases, r_beta = 0)
  d <- suppressMessages(runEpidemic(x$epop, x$sp))@dynamics
  i <- which(d$indCases == 1L)
  expect_true(all(d$Tinf[i] <= d$Tdeath[i]))
  
  # SEIR and SIDR: infection -> middle transition -> removal
  for (model in c("SEIR", "SIDR")) {
    x <- make_test_epidemic(model, group_list = group_list,
                            ind_cases = ind_cases, r_beta = 0)
    d <- suppressMessages(runEpidemic(x$epop, x$sp))@dynamics
    i <- which(d$indCases == 1L)
    
    expect_true(all(d$Tinf[i] <= d$Tsign[i]))
    expect_true(all(d$Tsign[i] <= d$Tdeath[i]))
  }
  
  # SEIDR: infection -> incubation -> detection -> removal
  x <- make_test_epidemic("SEIDR", group_list = group_list,
                          ind_cases = ind_cases, r_beta = 0)
  d <- suppressMessages(runEpidemic(x$epop, x$sp))@dynamics
  i <- which(d$indCases == 1L)
  
  expect_true(all(d$Tinf[i] <= d$Tinc[i]))
  expect_true(all(d$Tinc[i] <= d$Tsign[i]))
  expect_true(all(d$Tsign[i] <= d$Tdeath[i]))
})

test_that("runEpidemic validates its public inputs", {
  x <- make_test_epidemic("SIR")
  
  expect_warning(
    expect_error(
      runEpidemic(AlphaSimR::newPop(
        AlphaSimR::quickHaplo(nInd = 3, nChr = 1, segSites = 20),
        simParam = x$sp
      ), simParam = x$sp),
      "epop must be a PopEpidemic object"
    )
  )
  
  expect_error(
    runEpidemic(x$epop, simParam = AlphaSimR::SimParam$new(
      AlphaSimR::quickHaplo(nInd = 3, nChr = 1, segSites = 20)
    )),
    "simParam must be a SimParamEpidemic object"
  )
})

