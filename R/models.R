#' @title runEpidemic
#' @description
#' General wrapper for all compartmental models
#' @param epop A population with an epidemic, an object of class
#'  \code{\link{PopEpidemic-class}}
#' @param simParam a \code{\link{SimParamEpidemic}} object
#' @param ... Other inputs passed to specific functions. See details
#' @name runEpidemic
#' @rdname runEpidemic
#'
#' @examples
#' # Create founder haplotypes
#' founderPop <- quickHaplo(nInd = 10, nChr = 1, segSites = 10)
#'
#' # Set simulation parameters
#' SP <- SimParamEpidemic$new(founderPop, model = "SIR")
#' means <- rep(0, nchar(SP$model))
#' vars <- rep(1, nchar(SP$model))
#' SP$addTraitA(10, mean = means, var = vars, name = c("sus","inf","tol"))
#'
#' # Create population assuming one random index case
#' pop <- newPopEpidemic(founderPop)
#'
#' # set phenotypes
#' pop <- setPheno(pop, h2 = c(0.3, 0.4, 0.2))
#'
#' # Run Epidemic
#' pop <- runEpidemic(pop)
#'
#' @export
runEpidemic <- function(epop,simParam=NULL){
  
  if (!is(epop, "PopEpidemic") & is(epop, "Pop")){
    warning("epop is a Pop object. Use the function asPopEpidemic to convert it.")
  }
  
  if(is.null(simParam)) simParam = get("SP",envir=.GlobalEnv)
  
  stopifnot("simParam must be a SimParamEpidemic object"=
              is(simParam, "SimParamEpidemic"),
            "epop must be a PopEpidemic object."=
              is(epop, "PopEpidemic"))
  
  f <- get(paste0("model",toupper(simParam$model)))
  f(epop, simParam)
}

#' @title modelSIR
#' @description
#' Compartmental model SIR - Susceptibility, Infectivity, Recovrability
#'  (or Removed).
#'
#' @param epop A population with an epidemic, an object of class
#' \code{\link{PopEpidemic-class}}
#' @param simParam a \code{\link{SimParamEpidemic}} object
#'
#' @return Returns the same input object of type \code{\link{PopEpidemic-class}}
#' with the \code{dynamics} field populated with evolution of the epidemic
modelSIR <- function(epop, simParam) {
  
  et <- simParam$epi_traits
  
  #This function uses a cte val -> rempve and replace with a more efficient ver.
  .generate_sir_path <- function(epi_time, X, id, simParam) {
    simParam$removal_period * X[[et['t']]][[id]] + epi_time
# for a variable removal_period, gamma distribution is needed - deprecated
#rgamma(1L, simParam$RP_shape,
#       scale = simParam$RP_scale * X[[et['t']]][[id]]) + epi_time
}
  
  # Find the no. of individuals capable of infecting susceptible individuals at
  # some point (so includes those currently exposed and not yet infectious).
  .get_sir_infectives <- function(X) {
    X[, sum(status == "I")]
  }
  
  X <- cbind(data.table::as.data.table(epop@pheno),
             data.table::copy(epop@dynamics))
  X[, iid := epop@iid]
  
  purrr::walk(X[, .I[indCases == 1L]], \(i) {
    data.table::set(X, i, c("status", simParam$timings),
                    c(simParam$compartments[[2]], 0.0,
                      as.list(.generate_sir_path(0.0, X, i, simParam))))
  })
  
  #TODO: GE was Removed. Add back
  #
  Xgroups <- X |> split(by = "group")
  
  # Each group is independent, so it seems to be about ~33% faster to split
  # them by group, run them separately, then recombine them.
  
  Y <- purrr::map(Xgroups, \(X) {
    # This is a priority queue for the next event
    ni_events <- X[, .(.I, Tdeath)] |>
      data.table::melt(id.vars = "I", variable.name = "event",
                       value.name = "time") |>
      data.table::setorder(time, na.last = TRUE)
    ni_events[, event := NULL]
    
    # Start epidemic simulation loop ----
    epi_time <- 0.0
    while (.get_sir_infectives(X) > 0L) {
      
      # Calculate infection rates in each group ----
      #X[, group_inf := r_beta * GE * mean(inf * (status == "I"))]
      X[, group_inf := simParam$r_beta * mean(inf * (status == "I"))]
      
      # if S, infection at rate beta SI
      X[, inf_rate := sus * group_inf * (status == "S")]
      
      # id and time of next non-infection event
      t_next_event  <- ni_events$time[[1]]
      id_next_event <- ni_events$I[[1]]
      
      if (is.na(t_next_event)) t_next_event <- Inf
      
      # generate random timestep ----
      total_inf_rate <- sum(X$inf_rate)
      
      # calculate dt if infections event rate > 0
      dt <- rexp(1L) / total_inf_rate
      
      # check if next event is infection or non-infection ----
      if (epi_time + dt < t_next_event) {
        
        epi_time <- epi_time + dt
        
        # randomly select individual
        id_next_event <- sample(nrow(X),size = 1L,prob = X$inf_rate)
        
        xi <- X[, sample(x = .N, size = 1L, prob = inf * (status == "I"))]
        infd_by <- X$iid[[xi]]
        next_gen <- X$generation[[xi]] + 1L
        
        # get the recovery (death) time of the infected individual
        sir_path <- .generate_sir_path(epi_time, X, id_next_event, simParam)
        
        # store everything for that individual
        data.table::set(
          X, id_next_event,
          c("status", "Tinf", "Tdeath", "generation", "infected_by"),
          c("I", epi_time, as.list(sir_path), next_gen, infd_by))
        
        # also add the T_death to the priority queue
        ni_events[I == id_next_event, time := sir_path]
        
      } else {
        
        epi_time <- t_next_event
        
        # remove the death time for that (1st) individual from the queue
        data.table::set(ni_events, 1L, "time", NA)
        
        status <- X$status[[id_next_event]]
        
        if (status == "I") {
          data.table::set(X, id_next_event, "status", "R")
        } else {
          message("status = ", status)
          print(X[, .(group, indCases, status, Tinf, Tdeath,
                      group_inf, inf_rate)])
          print(X[id_next_event, .(group, indCases, status, Tinf, Tdeath,
                                   group_inf, inf_rate)])
          stop("selected ID ", id_next_event, "... unexpected event!")
          break
        }
      }
      
      data.table::setorder(ni_events, time, na.last = TRUE)
      
    }
    X
  }) |> rbindlist()
  
  # Restore epop's individual order. iid is NOT guaranteed to be sorted (it can
  # be arbitrary after crossing/subsetting), so we reorder to match epop@iid
  # rather than sorting ascending.
  Y <- Y[match(epop@iid, iid)]
  
  # remove iid though
  Y[, iid := NULL]
  
  # and the phenotypes
  Y[,colnames(epop@pheno) := NULL]
  
  # also these two columns are no longer needed
  Y[, c("group_inf", "inf_rate") := NULL]
  
  final_t <- max(Y$Tdeath, na.rm = TRUE)
  message(paste("- Final t =",final_t,", values are:",
                paste(capture.output(table(Y$status)), collapse = ", "), "\n"))
  epop@dynamics <- Y[]
  epop
  }

#' @title modelSEIR
#' @description
#' Compartmental model SEIR - Susceptible, Exposed, Infectious, Recovered
#'  (or Removed).
#'
#' @param epop A population with an epidemic, an object of class
#' \code{\link{PopEpidemic-class}}
#' @param simParam a \code{\link{SimParamEpidemic}} object
#'
#' @return Returns the same input object of type \code{\link{PopEpidemic-class}}
#'  with the \code{dynamics} field populated with evolution of the epidemic
modelSEIR <- function(epop, simParam) {
  
  et <- simParam$epi_traits
  
  # A susceptible individual's future trajectory is fixed at the moment of
  # exposure. Returns the cumulative event times [Tsign, Tdeath] measured from
  # `epi_time`: first the E -> I transition (latent period, scaled by the `lat`
  # phenotype), then the I -> R transition (infectious period, scaled by `tol`).
  # this function's now using constant values -> TODO: remove and replace
  .generate_seir_path <- function(epi_time, X, id, simParam) {
    #rgamma(2L, c(simParam$LP_shape, simParam$RP_shape),
    #       scale = c(simParam$LP_scale * X[[et['l']]][[id]],
    #                 simParam$RP_scale * X[[et['t']]][[id]])) |> cumsum() +
    #  epi_time
    c(simParam$latent_period, simParam$removal_period) * 
      c(X[[et['l']]][[id]], X[[et['t']]][[id]]) |> cumsum() + epi_time
  }
  
  # Individuals that can still drive the epidemic, i.e. those that are
  # infectious ("I") or will become infectious ("E"). Once none remain the
  # epidemic is over.
  .get_seir_infectives <- function(X) {
    X[, sum(status %in% c("E", "I"))]
  }
  
  X <- cbind(data.table::as.data.table(epop@pheno),
             data.table::copy(epop@dynamics))
  X[, iid := epop@iid]
  
  # Seed indCases into the E compartment with a fixed future trajectory
  purrr::walk(X[, .I[indCases == 1L]], \(i) {
    data.table::set(X, i, c("status", simParam$timings),
                    c(simParam$compartments[[2]], 0.0,
                      as.list(.generate_seir_path(0.0, X, i, simParam))))
  })
  
  #TODO: GE was Removed. Add back
  
  Xgroups <- X |> split(by = "group")
  
  # Each group is independent, so it seems to be about ~33% faster to split
  # them by group, run them separately, then recombine them.
  
  Y <- purrr::map(Xgroups, \(X) {
    # Priority queue for the next non-infection event. Unlike SIR, each
    # individual now has TWO scheduled non-infection events (E -> I at Tsign and
    # I -> R at Tdeath), so the melt produces two rows per individual. The
    # transition that is actually applied is decided from the individual's
    # status when its event is popped.
    ni_events <- X[, .(.I, Tsign, Tdeath)] |>
      data.table::melt(id.vars = "I", variable.name = "event",
                       value.name = "time") |>
      data.table::setorder(time, na.last = TRUE)
    ni_events[, event := NULL]
    
    # Start epidemic simulation loop ----
    epi_time <- 0.0
    while (.get_seir_infectives(X) > 0L) {
      
      # Calculate infection rates in each group ----
      # Only "I" individuals are infectious; "E" individuals are not yet.
      #X[, group_inf := r_beta * GE * mean(inf * (status == "I"))]
      X[, group_inf := simParam$r_beta * mean(inf * (status == "I"))]
      
      # if S, infection at rate beta SI
      X[, inf_rate := sus * group_inf * (status == "S")]
      
      # id and time of next non-infection event
      t_next_event  <- ni_events$time[[1]]
      id_next_event <- ni_events$I[[1]]
      
      if (is.na(t_next_event)) t_next_event <- Inf
      
      # generate random timestep ----
      total_inf_rate <- sum(X$inf_rate)
      
      # calculate dt if infections event rate > 0
      dt <- rexp(1L) / total_inf_rate
      
      # check if next event is infection or non-infection ----
      if (epi_time + dt < t_next_event) {
        
        epi_time <- epi_time + dt
        
        # randomly select the individual that becomes infected
        id_next_event <- sample(nrow(X), size = 1L, prob = X$inf_rate)
        
        # randomly select the infector (only "I" individuals are infectious)
        xi <- X[, sample(x = .N, size = 1L, prob = inf * (status == "I"))]
        infd_by <- X$iid[[xi]]
        next_gen <- X$generation[[xi]] + 1L
        
        seir_path <- .generate_seir_path(epi_time, X, id_next_event, simParam)
        
        data.table::set(
          X, id_next_event,
          c("status", "Tinf", "Tsign", "Tdeath", "generation", "infected_by"),
          c("E", epi_time, as.list(seir_path), next_gen, infd_by))
        
        # Update this individual's two queued events. seir_path has length 2 and
        # matches the two rows for this individual; the assignment order is
        # irrelevant since both times belong to the same individual and are
        # resolved by status when popped.
        ni_events[I == id_next_event, time := seir_path]
        
      } else {
        
        epi_time <- t_next_event
        
        data.table::set(ni_events, 1L, "time", NA)
        
        status <- X$status[[id_next_event]]
        
        if (status == "E") {
          data.table::set(X, id_next_event, "status", "I")
        } else if (status == "I") {
          data.table::set(X, id_next_event, "status", "R")
        } else {
          message("status = ", status)
          print(X[, .(group, indCases, status, Tinf, Tsign, Tdeath,
                      group_inf, inf_rate)])
          print(X[id_next_event, .(group, indCases, status, Tinf, Tsign, Tdeath,
                                   group_inf, inf_rate)])
          stop("selected ID ", id_next_event, "... unexpected event!")
          break
        }
      }
      
      data.table::setorder(ni_events, time, na.last = TRUE)
      
    }
    X
  }) |> rbindlist()
  
  # Restore epop's individual order. iid is NOT guaranteed to be sorted (it can
  # be arbitrary after crossing/subsetting), so we reorder to match epop@iid
  # rather than sorting ascending.
  Y <- Y[match(epop@iid, iid)]
  
  # remove iid though
  Y[, iid := NULL]
  
  # and the phenotypes
  Y[,colnames(epop@pheno) := NULL]
  
  # also these two columns are no longer needed
  Y[, c("group_inf", "inf_rate") := NULL]
  
  final_t <- max(Y$Tdeath, na.rm = TRUE)
  message(paste("- Final t =",final_t,", values are:",
                paste(capture.output(table(Y$status)), collapse = ", "), "\n"))
  
  epop@dynamics <- Y[]
  epop
}

#' @title modelSIDR
#' @description
#' Compartmental model SIDR - Susceptible, Infectious, Detected, Recovered
#'  (or Removed).
#'
#' @param epop A population with an epidemic, an object of class
#' \code{\link{PopEpidemic-class}}
#' @param simParam a \code{\link{SimParamEpidemic}} object
#'
#' @return Returns the same input object of type \code{\link{PopEpidemic-class}}
#'  with the \code{dynamics} field populated with evolution of the epidemic
modelSIDR <- function(epop, simParam) {
  
  et <- simParam$epi_traits
  
  # A susceptible individual's future trajectory is fixed at the moment of
  # infection. Returns the cumulative event times [Tsign, Tdeath] measured from
  # `epi_time`: first the I -> D transition (detection period, scaled by the
  # `det` phenotype), then the D -> R transition (removal period, scaled by
  # `tol`).
  # This function now uses cte values. TODO: remove and replace
  .generate_sidr_path <- function(epi_time, X, id, simParam) {
    #cumsum(rgamma(2L, c(simParam$DP_shape, simParam$RP_shape),
    #              scale = c(simParam$DP_scale * X[[et['d']]][[id]],
    #                        simParam$RP_scale * X[[et['t']]][[id]]))) +
    #  epi_time
    c(simParam$detection_period, simParam$removal_period) *
      c(X[[et['d']]][[id]], X[[et['t']]][[id]]) |> cumsum() + epi_time
  }
  
  # Individuals that can still drive the epidemic. In SIDR BOTH "I" and "D" are
  # infectious. Once none remain the epidemic is over.
  .get_sidr_infectives <- function(X) {
    X[, sum(status %in% c("I", "D"))]
  }
  
  X <- cbind(data.table::as.data.table(epop@pheno),
             data.table::copy(epop@dynamics))
  X[, iid := epop@iid]
  
  # Seed indCases into the I compartment with a fixed future trajectory
  purrr::walk(X[, .I[indCases == 1L]], \(i) {
    data.table::set(X, i, c("status", simParam$timings),
                    c(simParam$compartments[[2]], 0.0,
                      as.list(.generate_sidr_path(0.0, X, i, simParam))))
  })
  
  #TODO: GE was Removed. Add back
  
  Xgroups <- X |> split(by = "group")
  
  # Each group is independent, so it seems to be about ~33% faster to split
  # them by group, run them separately, then recombine them.
  
  Y <- purrr::map(Xgroups, \(X) {
    # Priority queue for the next non-infection event. As with SEIR each
    # individual has TWO scheduled non-infection events (I -> D at Tsign and
    # D -> R at Tdeath), so the melt produces two rows per individual. The
    # transition that is actually applied is decided from the individual's
    # status when its event is popped.
    ni_events <- X[, .(.I, Tsign, Tdeath)] |>
      data.table::melt(id.vars = "I", variable.name = "event",
                       value.name = "time") |>
      data.table::setorder(time, na.last = TRUE)
    ni_events[, event := NULL]
    
    # Start epidemic simulation loop ----
    epi_time <- 0.0
    while (.get_sidr_infectives(X) > 0L) {
      # Calculate infection rates in each group ----
      # Both "I" and "D" individuals are infectious.
      #X[, group_inf := r_beta * GE * mean(inf * (status %in% c("I", "D")))]
      X[, group_inf := simParam$r_beta * mean(inf * (status %in% c("I", "D")))]
      
      # if S, infection at rate beta SI
      X[, inf_rate := sus * group_inf * (status == "S")]
      
      # id and time of next non-infection event
      t_next_event  <- ni_events$time[[1]]
      id_next_event <- ni_events$I[[1]]
      
      if (is.na(t_next_event)) t_next_event <- Inf
      
      # generate random timestep ----
      total_inf_rate <- sum(X$inf_rate)
      
      # calculate dt if infections event rate > 0
      dt <- rexp(1L) / total_inf_rate
      
      # check if next event is infection or non-infection ----
      if (epi_time + dt < t_next_event) {
        
        epi_time <- epi_time + dt
        
        # randomly select the individual that becomes infected
        id_next_event <- sample(nrow(X), size = 1L, prob = X$inf_rate)
        
        # randomly select the infector (both "I" and "D" are infectious)
        xi <- X[, sample(x = .N, size = 1L,
                         prob = inf * (status %in% c("I", "D")))]
        infd_by <- X$iid[[xi]]
        next_gen <- X$generation[[xi]] + 1L
        
        sidr_path <- .generate_sidr_path(epi_time, X, id_next_event, simParam)
        
        set(X, id_next_event,
            c("status", "Tinf", "Tsign", "Tdeath", "generation", "infected_by"),
            c("I", epi_time, as.list(sidr_path), next_gen, infd_by))
        
        # Update this individual's two queued events.
        ni_events[I == id_next_event, time := sidr_path]
        
      } else {
        
        epi_time <- t_next_event
        
        data.table::set(ni_events, 1L, "time", NA)
        
        status <- X$status[[id_next_event]]
        
        if (status == "I") {
          data.table::set(X, id_next_event, "status", "D")
        } else if (status == "D") {
          data.table::set(X, id_next_event, "status", "R")
        } else {
          message("status = ", status)
          print(X[, .(group, indCases, status, Tinf, Tsign, Tdeath,
                      group_inf, inf_rate)])
          print(X[id_next_event, .(group, indCases, status, Tinf, Tsign, Tdeath,
                                   group_inf, inf_rate)])
          stop("selected ID ", id_next_event, "... unexpected event!")
          break
        }
      }
      
      data.table::setorder(ni_events, time, na.last = TRUE)
      
    }
    X
  }) |> rbindlist()
  
  # Restore epop's individual order.
  Y <- Y[match(epop@iid, iid)]
  
  # remove iid though
  Y[, iid := NULL]
  
  # and the phenotypes
  Y[,colnames(epop@pheno) := NULL]
  
  # also these two columns are no longer needed
  Y[, c("group_inf", "inf_rate") := NULL]
  
  final_t <- max(Y$Tdeath, na.rm = TRUE)
  message(sprintf("- Final t = %f, values are:", final_t),
          paste(capture.output(table(Y$status)), collapse = ", "), "\n")
  
  epop@dynamics <- Y[]
  epop
}


#' @title modelSEIDR
#' @description
#' Compartmental model SEIDR - Susceptible, Exposed, Infectious, Detected,
#'  Recovered (or Removed).
#'
#' @param epop A population with an epidemic, an object of class
#' \code{\link{PopEpidemic-class}}
#' @param simParam a \code{\link{SimParamEpidemic}} object
#'
#' @return Returns the same input object of type \code{\link{PopEpidemic-class}}
#'  with the \code{dynamics} field populated with evolution of the epidemic
modelSEIDR <- function(epop, simParam) {
  
  et <- simParam$epi_traits
  
  # A susceptible individual's future trajectory is fixed at the moment of
  # exposure. Returns the cumulative event times [Tinc, Tsign, Tdeath] measured
  # from `epi_time`: E -> I (latent period, scaled by `lat`), then I -> D
  # (detection period, scaled by `det`), then D -> R (removal period, scaled by
  # `tol`).
  #This function now uses cte values -> todo: remove and replace
  .generate_seidr_path <- function(epi_time, X, id, simParam) {
    #cumsum(rgamma(3L,
    #              c(simParam$LP_shape, simParam$DP_shape, simParam$RP_shape),
    #              scale = c(simParam$LP_scale * X[[et['l']]][[id]],
    #                        simParam$DP_scale * X[[et['d']]][[id]],
    #                        simParam$RP_scale * X[[et['t']]][[id]]))) + epi_time
    c(simParam$latent_period,simParam$detection_period,simParam$removal_period)*
      c(X[[et['l']]][[id]], X[[et['d']]][[id]],X[[et['t']]][[id]]) |> 
      cumsum() + epi_time
  }
  
  # Individuals that can still drive the epidemic, i.e. those that are
  # infectious ("I"/"D") or will become infectious ("E"). Once none remain the
  # epidemic is over.
  .get_seidr_infectives <- function(X) {
    X[, sum(status %in% c("E", "I", "D"))]
  }
  
  X <- cbind(data.table::as.data.table(epop@pheno),
             data.table::copy(epop@dynamics))
  X[, iid := epop@iid]
  
  # Seed indCases into the E compartment with a fixed future trajectory
  purrr::walk(X[, .I[indCases == 1L]], \(i) {
    data.table::set(X, i, c("status", simParam$timings),
                    c("E", 0.0,
                      as.list(.generate_seidr_path(0.0, X, i, simParam))))
  })
  
  #TODO: GE was Removed. Add back
  
  Xgroups <- X |> split(by = "group")
  
  # Each group is independent, so it seems to be about ~33% faster to split
  # them by group, run them separately, then recombine them.
  
  Y <- purrr::map(Xgroups, \(X) {
    # Priority queue for the next non-infection event. Each individual now has
    # THREE scheduled non-infection events (E -> I at Tinc, I -> D at Tsign,
    # D -> R at Tdeath), so the melt produces three rows per individual. The
    # transition that is actually applied is decided from the individual's
    # status when its event is popped.
    ni_events <- X[, .(.I, Tinc, Tsign, Tdeath)] |>
      data.table::melt(id.vars = "I", variable.name = "event",
                       value.name = "time") |>
      data.table::setorder(time, na.last = TRUE)
    ni_events[, event := NULL]
    
    # Start epidemic simulation loop ----
    epi_time <- 0.0
    while (.get_seidr_infectives(X) > 0L) {
      
      # Calculate infection rates in each group ----
      # Both "I" and "D" individuals are infectious; "E" individuals are not.
      #X[, group_inf := r_beta * GE * mean(inf * (status %in% c("I", "D")))]
      X[, group_inf := simParam$r_beta * mean(inf * (status %in% c("I", "D")))]
      
      # if S, infection at rate beta SI
      X[, inf_rate := sus * group_inf * (status == "S")]
      
      # id and time of next non-infection event
      t_next_event  <- ni_events$time[[1]]
      id_next_event <- ni_events$I[[1]]
      
      if (is.na(t_next_event)) t_next_event <- Inf
      
      # generate random timestep ----
      total_inf_rate <- sum(X$inf_rate)
      
      # calculate dt if infections event rate > 0
      dt <- rexp(1L) / total_inf_rate
      
      # check if next event is infection or non-infection ----
      if (epi_time + dt < t_next_event) {
        
        epi_time <- epi_time + dt
        
        # randomly select the individual that becomes infected
        id_next_event <- sample(nrow(X), size = 1L, prob = X$inf_rate)
        
        # randomly select the infector (both "I" and "D" are infectious)
        xi <- X[, sample(x = .N, size = 1L,
                         prob = inf * (status %in% c("I", "D")))]
        infd_by <- X$iid[[xi]]
        next_gen <- X$generation[[xi]] + 1L
        
        seidr_path <- .generate_seidr_path(epi_time, X, id_next_event, simParam)
        
        data.table::set(
          X, id_next_event,
          c("status", "Tinf", "Tinc", "Tsign", "Tdeath", "generation",
            "infected_by"),
          c("E", epi_time, as.list(seidr_path), next_gen, infd_by))
        
        # Update this individual's three queued events. seidr_path has length 3
        # and matches the three rows for this individual; the assignment order
        # is irrelevant since all times belong to the same individual and are
        # resolved by status when popped.
        ni_events[I == id_next_event, time := seidr_path]
        
      } else {
        
        epi_time <- t_next_event
        
        data.table::set(ni_events, 1L, "time", NA)
        
        status <- X$status[[id_next_event]]
        
        if (status == "E") {
          data.table::set(X, id_next_event, "status", "I")
        } else if (status == "I") {
          data.table::set(X, id_next_event, "status", "D")
        } else if (status == "D") {
          data.table::set(X, id_next_event, "status", "R")
        } else {
          message("status = ", status)
          print(X[, .(group, indCases, status, Tinf, Tinc, Tsign, Tdeath,
                      group_inf, inf_rate)])
          print(X[id_next_event, .(group, indCases, status, Tinf, Tinc, Tsign,
                                   Tdeath, group_inf, inf_rate)])
          stop("selected ID ", id_next_event, "... unexpected event!")
          break
        }
      }
      
      data.table::setorder(ni_events, time, na.last = TRUE)
      
    }
    X
  }) |> rbindlist()
  
  # Restore epop's individual order.
  Y <- Y[match(epop@iid, iid)]
  
  # remove iid though
  Y[, iid := NULL]
  
  # and the phenotypes
  Y[,colnames(epop@pheno) := NULL]
  
  # also these two columns are no longer needed
  Y[, c("group_inf", "inf_rate") := NULL]
  
  final_t <- max(Y$Tdeath, na.rm = TRUE)
  message(sprintf("- Final t = %f, values are:", final_t),
          paste(capture.output(table(Y$status)), collapse = ", "), "\n")
  
  epop@dynamics <- Y[]
  epop
}


#' @title modelSIS
#' @description
#' Compartmental model SIS - Susceptible, Infectious, Susceptible again
#'
#' @param epop A population with an epidemic, an object of class
#' \code{\link{PopEpidemic-class}}
#' @param simParam a \code{\link{SimParamEpidemic}} object
#'
#' @return Returns the same input object of type \code{\link{PopEpidemic-class}}
#'  with the \code{dynamics} field populated with evolution of the epidemic
modelSIS <- function(epop, simParam) {
  
  et <- simParam$epi_traits
  
  # A Susceptible individual's future disease trajectory is fixed at the point of
  # exposure.
  generate_sis_path <- function(epi_time, X, id, simParam) {
    simParam$removal_period * X[[et['t']]][[id]] + epi_time
    #for a variable removal_period, gamma distribution is needed - deprecated
    #rgamma(1L, simParam$RP_shape,
      #scale = simParam$RP_scale * X[[et['t']]][[id]]) + epi_time
  }
  
  r_beta <- simParam$r_beta
  t_end  <- simParam$t_final

  X <- cbind(data.table::as.data.table(epop@pheno),
             data.table::copy(epop@dynamics))

  nax <- function(x, y) {
    if (length(x) == 1 && is.na(x[[1]])) y else c(x, y)
  }

  X[, iid := epop@iid]

  purrr::walk(X[, .I[indCases == 1L]], \(i) {
    data.table::set(X, i, c("status", simParam$timings),
                    c(simParam$compartments[[2]], 0.0,
                      as.list(.generate_sir_path(0.0, X, i, simParam))))
  })
  
  Xgroups <- X |> split(by = "group")
  
  Y <- map(Xgroups, \(X) {
    
    # This is unique to the SIS and SIRS models
    X[, `:=`(Tinf = as.list(Tinf),
             Tdeath = as.list(Tdeath),
             generation = as.list(generation),
             infected_by = as.list(infected_by))]
    
    # This is a priority queue for the next event
    ni_events <- X[, .(.I, map_dbl(Tdeath, last))] |>
      data.table::melt(id.vars = "I", variable.name = "event", 
      value.name = "time") |>
      data.table::setorder(time, na.last = TRUE)
    ni_events[, event := NULL]
    
    # Start epidemic simulation loop
    epi_time <- 0.0
    while (epi_time < t_end && X[, sum(status == "I") > 0L]) {
      
      # Calculate infection rates in each group
      # this is the sum of the log infectivitities
      #X[, group_inf := r_beta * GE * mean(inf * (status == "I"))]
      X[, group_inf := r_beta * mean(inf * (status == "I"))]
      
      # if S, infection at rate beta SI
      X[, inf_rate := sus * group_inf * (status == "S")]
      
      # id and time of next non-infection event
      t_next_event  <- ni_events$time[[1]]
      id_next_event <- ni_events$I[[1]]
      
      # generate random timestep -
      total_inf_rate <- sum(X$inf_rate)
      
      # calculate dt if infections event rate > 0
      dt <- rexp(1L) / total_inf_rate
      
      # check if next event is infection or non-infection
      if (epi_time + dt < t_next_event) {
        
        epi_time <- epi_time + dt
        
        # randomly select individual
        id_next_event <- sample(nrow(X),
                                size = 1L,
                                prob = X$inf_rate)
        
        xi <- X[, sample(.N, 1L, prob = inf * (status == "I"))]
        infd_by <- X$iid[[xi]]
        next_gen <- X$generation[[xi]] + 1L
        
        sis_path <- generate_sis_path(epi_time, X, id_next_event, params)
        
        set(X, id_next_event,
            c("status", "Tinf", "Tdeath", "generation", "infected_by"),
            list("I",
                 nax(X$Tinf[[id_next_event]], epi_time),
                 nax(X$Tdeath[[id_next_event]], sis_path),
                 nax(X$generation[[id_next_event]], next_gen),
                 nax(X$infected_by[[id_next_event]], infd_by)))
        
        ni_events[I == id_next_event, time := sis_path]
        
      } else {
        
        epi_time <- t_next_event
        
        data.table::set(ni_events, 1L, "time", NA)
        
        status <- X[id_next_event, status]
        
        if (status == "I") {
          data.table::set(X, id_next_event, "status", "S")
        } else {
          print(X[, .(group, donor, status, Tinf, Tdeath, group_inf, inf_rate)])
          print(X[id_next_event, .(group, donor, status, Tinf, Tdeath, group_inf, inf_rate)])
          stop(str_c("selected ID ", id_next_event, "... unexpected event!"))
          break
        }
      }
      
      # Every event we set the order of the ni_events. This should be fast
      # though as data.table uses a sensible sort method.
      data.table::setorder(ni_events, time, na.last = TRUE)
    }
    X
  }) |> rbindlist()
  
  final_t <- Y$Tdeath |> list_c() |> max(na.rm = TRUE) |> signif(5)
  
  message(str_glue("- Final t = {final_t}, values are:"),
          str_flatten(capture.output(table(Y$status)), "\n"))
  
  Y <- Y[match(epop@iid, iid)]

  # remove iid though
  Y[, iid := NULL]
  
  # and the phenotypes
  Y[,colnames(epop@pheno) := NULL]
  
  # also these two columns are no longer needed
  Y[, c("group_inf", "inf_rate") := NULL]
  
  epop@dynamics <- Y[]
  
  return(epop)
}