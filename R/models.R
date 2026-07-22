#' @title run
#' @description
#' General wrapper for all compartmental model
#' @param epop A population with an epidemy, an object of class
#' \code{\link{PopEpidemic-class}}
#' @param simParam a \code{\link{SimParamEpidemic}} object
#' @name run
#' @rdname run
#' @export
run <- function(epop,simParam=NULL){
  if(is.null(simParam)) simParam = get("SP",envir=.GlobalEnv)
  f <- get(paste0("model",toupper(simParam$model)))
  f(epop, simParam)
}

#' @title modelSIR
#' @description
#' Compartmental model SIR - Susceptibility, Infectivity, Recovrability
#'  (or Removed)).
#'
#' @param epop A population with an epidemy, an object of class
#' \code{\link{PopEpidemic-class}}
#' @param simParam a \code{\link{SimParamEpidemic}} object
#'
#' @return Returns the same input object of type \code{\link{PopEpidemic-class}} with
#'  the \code{dynamic} field populated with evolution of the epidemy
modelSIR <- function(epop, simParam = NULL) {
  DEBUG <- F
  if(is.null(simParam)) simParam = get("SP",envir=.GlobalEnv)

  stopifnot("simParam must be a SimParamEpidemic object"=
              is(simParam, "SimParamEpidemic"),
            "epop must be a PopEpidemic object"=
              is(epop, "PopEpidemic"),
            "Not an SIR model"=
              toupper(simParam$model)=="SIR")

  .generate_sir_path <- function(epi_time, X, id, simParam) {
    Tinf   <- epi_time
    Tdeath <- Tinf + rgamma(1L, simParam$RP_shape,
                            scale = simParam$RP_scale * X$tol[[id]])
    list("I", Tinf, Tdeath)
  }

  # Find the no. of individuals capable of infecting susceptible individuals at
  # some point (so includes those currently exposed and not yet infectious).
  .get_sir_infectives <- function(X) {
    X[, sum(status == "I")]
  }

  # Find the time of the next non-infection event, and the id of the individual
  .next_sir_ni_event <- function(X, epi_time) {
    as.list(X[, .(.I,
                  Tmin = fifelse(Tdeath > epi_time, Tdeath, Inf))]
            [, list(t_next_event = min(Tmin, na.rm = TRUE),
                    id_next_event = which.min(Tmin))])
  }

  X <- .epi_working_table(epop)
  .epi_init_donors(X, .generate_sir_path, simParam)

  #TODO: GE was Removed. Add back
  #
  Xgroups <- X |> split(by = "group")

  # Each group is independent, so it seems to be about ~33% faster to split
  # them by group, run them separately, then recombine them.

  Y <- map(Xgroups, \(X) {
    # Start epidemic simulation loop ----
    epi_time <- 0.0
    while (.get_sir_infectives(X) > 0L) {
      if (DEBUG) message("time = ", signif(epi_time, 5))

      # Calculate infection rates in each group ----
      #X[, group_inf := r_beta * GE * mean(fifelse(status == "I", inf, 0.0))]
      X[, group_inf := simParam$r_beta * mean(fifelse(status == "I", inf, 0.0))]

      # if S, infection at rate beta SI
      X[, event_rate := fifelse(status == "S", sus * group_inf, 0.0)]

      # id and time of next non-infection event
      ni_event      <- .next_sir_ni_event(X, epi_time)
      t_next_event  <- ni_event$t_next_event
      id_next_event <- ni_event$id_next_event

      if (is.na(t_next_event)) t_next_event <- Inf

      if (DEBUG) message("Next NI event id = ", id_next_event, " at t = ",
                         signif(t_next_event, 5))

      # generate random timestep ----
      total_event_rate <- sum(X$event_rate)

      # calculate dt if infections event rate > 0
      dt <- if (total_event_rate > 0.0) {
        rexp(1L, rate = total_event_rate)
      } else {
        Inf
      }

      if (DEBUG) message("Total infections event rate = ",
                         signif(total_event_rate, 5))

      # check if next event is infection or non-infection ----
      if (epi_time + dt < t_next_event) {
        if (DEBUG) message("next event is infection at t = ", signif(epi_time, 5))

        epi_time <- epi_time + dt

        # randomly select individual
        id_next_event <- sample(nrow(X),
                                size = 1L,
                                prob = X$event_rate)

        infectives <- X[status == "I", .(.I, status, inf)]
        infd_by <- .safe_sample(x = infectives$I,
                                size = 1L,
                                prob = infectives$inf)
        next_gen <- X$generation[[infd_by]] + 1L

        set(X, id_next_event, c("status", "Tinf", "Tdeath"),
            .generate_sir_path(epi_time, X, id_next_event, simParam))
        set(X, id_next_event, c("generation", "infected_by"),
            list(next_gen, infd_by))

        if (DEBUG) message("ID ", id_next_event,
                           ": S -> I, infected by ID ", infd_by)
      } else {
        if (DEBUG) message("next event is non-infection at t = ",
                           signif(t_next_event, 5))

        epi_time <- t_next_event

        status <- X$status[[id_next_event]]

        if (status == "I") {
          if (DEBUG) message("ID ", id_next_event, ": I -> R")
          set(X, id_next_event, "status", "R")
        } else {
          message("status = ", status)
          print(X[, .(group, donor, status, Tinf, Tdeath,
                      group_inf, event_rate)])
          print(X[id_next_event, .(group, donor, status, Tinf, Tdeath,
                                   group_inf, event_rate)])
          stop("selected ID ", id_next_event, "... unexpected event!")
          break
        }
      }

      if (DEBUG == 2) {
        print(X[,.(group, donor, status, Tinf, Tdeath, group_inf, event_rate)])
      }
    }
    X
  }) |> rbindlist()

  .finalise_epi_dynamics(epop, Y)
}

#' @title modelSIDR
#' @description
#' Compartmental model SIDR - Susceptibility, Infectivity, Detectability,
#'  Recoverability (or Removed). Both the undetected (I) and detected (D)
#'  individuals are infectious.
#'
#' @param epop A population with an epidemy, an object of class
#' \code{\link{PopEpidemic-class}}
#' @param simParam a \code{\link{SimParamEpidemic}} object
#'
#' @return Returns the same input object of type \code{\link{PopEpidemic-class}}
#'  with the \code{dynamic} field populated with evolution of the epidemy
modelSIDR <- function(epop, simParam = NULL) {
  DEBUG <- F
  if (is.null(simParam)) simParam = get("SP", envir = .GlobalEnv)
  
  stopifnot("simParam must be a SimParamEpidemic object" =
              is(simParam, "SimParamEpidemic"),
            "epop must be a PopEpidemic object" =
              is(epop, "PopEpidemic"),
            "Not an SIDR model" =
              toupper(simParam$model) == "SIDR")
  
  # A susceptible individual's future trajectory is fixed at infection:
  # I (Tinf) -> D (Tsym) -> R (Tdeath). Detection uses 'det', removal uses 'tol'.
  .generate_sidr_path <- function(epi_time, X, id, simParam) {
    Tinf   <- epi_time
    Tsym   <- Tinf + rgamma(1L, simParam$DP_shape,
                            scale = simParam$DP_scale * X$det[[id]])
    Tdeath <- Tsym + rgamma(1L, simParam$RP_shape,
                            scale = simParam$RP_scale * X$tol[[id]])
    list("I", Tinf, Tsym, Tdeath)
  }
  
  # In SIDR both undetected (I) and detected (D) individuals are infectious.
  .get_sidr_infectives <- function(X) {
    X[, sum(status %in% c("I", "D"))]
  }
  
  # Next non-infection event. Which transition is pending depends on the current
  # status: I -> D fires at Tsym, D -> R fires at Tdeath. (epi_time kept for
  # signature parity with the SIR helper; status already disambiguates here.)
  .next_sidr_ni_event <- function(X, epi_time) {
    as.list(X[, .(.I,
                  Tmin = fcase(status == "I", Tsym,
                               status == "D", Tdeath,
                               default = Inf))]
            [, list(t_next_event = min(Tmin, na.rm = TRUE),
                    id_next_event = which.min(Tmin))])
  }
  
  #TODO: Quickest way to reuse Jamie's code but maybe not memory friendly
  X <- epop@pheno |> as.data.table() |> cbind(epop@dynamics)
  
  # Fix the disease trajectory for the initial donors
  purrr::walk(X[, .I[donor == 1L]], \(i) {
    data.table::set(X, i, c("status", simParam$timings),
                    .generate_sidr_path(0.0, X, i, simParam))
  })
  
  #TODO: GE was Removed. Add back
  Xgroups <- X |> split(by = "group")
  
  # Each group is independent, so it seems to be about ~33% faster to split
  # them by group, run them separately, then recombine them.
  Y <- map(Xgroups, \(X) {
    # Start epidemic simulation loop ----
    epi_time <- 0.0
    while (.get_sidr_infectives(X) > 0L) {
      if (DEBUG) message("time = ", signif(epi_time, 5))
      
      # Calculate infection rate in the group (both I and D contribute) ----
      X[, group_inf := simParam$r_beta *
          mean(fifelse(status %in% c("I", "D"), inf, 0.0))]
      
      # if S, infection at rate beta SI
      X[, event_rate := fifelse(status == "S", sus * group_inf, 0.0)]
      
      # id and time of next non-infection event
      ni_event      <- .next_sidr_ni_event(X, epi_time)
      t_next_event  <- ni_event$t_next_event
      id_next_event <- ni_event$id_next_event
      
      if (is.na(t_next_event)) t_next_event <- Inf
      
      if (DEBUG) message("Next NI event id = ", id_next_event, " at t = ",
                         signif(t_next_event, 5))
      
      # generate random timestep ----
      total_event_rate <- sum(X$event_rate)
      
      dt <- if (total_event_rate > 0.0) {
        rexp(1L, rate = total_event_rate)
      } else {
        Inf
      }
      
      if (DEBUG) message("Total infections event rate = ",
                         signif(total_event_rate, 5))
      
      # check if next event is infection or non-infection ----
      if (epi_time + dt < t_next_event) {
        if (DEBUG) message("next event is infection at t = ", signif(epi_time, 5))
        
        epi_time <- epi_time + dt
        
        # randomly select individual to be infected
        id_next_event <- sample(nrow(X), size = 1L, prob = X$event_rate)
        
        # select the infector, weighted by infectivity (from I or D)
        infectives <- X[status %in% c("I", "D"), .(.I, status, inf)]
        infd_by <- .safe_sample(x = infectives$I, size = 1L,
                                prob = infectives$inf)
        next_gen <- X$generation[[infd_by]] + 1L
        
        set(X, id_next_event, c("status", "Tinf", "Tsym", "Tdeath"),
            .generate_sidr_path(epi_time, X, id_next_event, simParam))
        set(X, id_next_event, c("generation", "infected_by"),
            list(next_gen, infd_by))
        
        if (DEBUG) message("ID ", id_next_event,
                           ": S -> I, infected by ID ", infd_by)
      } else {
        if (DEBUG) message("next event is non-infection at t = ",
                           signif(t_next_event, 5))
        
        epi_time <- t_next_event
        
        status <- X$status[[id_next_event]]
        
        if (status == "I") {
          if (DEBUG) message("ID ", id_next_event, ": I -> D")
          set(X, id_next_event, "status", "D")
        } else if (status == "D") {
          if (DEBUG) message("ID ", id_next_event, ": D -> R")
          set(X, id_next_event, "status", "R")
        } else {
          message("status = ", status)
          print(X[, .(group, donor, status, Tinf, Tsym, Tdeath,
                      group_inf, event_rate)])
          print(X[id_next_event, .(group, donor, status, Tinf, Tsym, Tdeath,
                                   group_inf, event_rate)])
          stop("selected ID ", id_next_event, "... unexpected event!")
          break
        }
      }
      
      if (DEBUG == 2) {
        print(X[, .(group, donor, status, Tinf, Tsym, Tdeath,
                    group_inf, event_rate)])
      }
    }
    X
  }) |> rbindlist()
  
  final_t <- max(Y$Tdeath, na.rm = TRUE) |> signif(5)
  
  message(sprintf("- Final t = %f, values are:", final_t),
          paste(capture.output(table(Y$status)), collapse = ", "), "\n")
  
  # tidy up X
  Y[, c("group_inf", "event_rate") := NULL]
  Y[, parasites := !is.na(Tinf)] # what is this column
  
  epop@dynamics <- Y
  
  epop
}

#' @title modelSEIR
#' @description
#' Compartmental model SEIR - Susceptibility, Exposed (latent), Infectivity,
#'  Recoverability (or Removed). Exposed (E) individuals are infected but NOT
#'  yet infectious; only infectious (I) individuals can transmit.
#'
#' @param epop A population with an epidemy, an object of class
#' \code{\link{PopEpidemic-class}}
#' @param simParam a \code{\link{SimParamEpidemic}} object
#'
#' @return Returns the same input object of type \code{\link{PopEpidemic-class}}
#'  with the \code{dynamic} field populated with evolution of the epidemy
modelSEIR <- function(epop, simParam = NULL) {
  DEBUG <- F
  if (is.null(simParam)) simParam = get("SP", envir = .GlobalEnv)
  
  stopifnot("simParam must be a SimParamEpidemic object" =
              is(simParam, "SimParamEpidemic"),
            "epop must be a PopEpidemic object" =
              is(epop, "PopEpidemic"),
            "Not an SEIR model" =
              toupper(simParam$model) == "SEIR")
  
  # A susceptible individual's future trajectory is fixed at infection:
  # E (Tinf) -> I (Tsym) -> R (Tdeath). The latent period E->I uses 'lat'
  # (LP_*); the removal period I->R uses 'tol' (RP_*). NB: for SEIR "Tsym" is
  # the onset-of-infectiousness time (end of latency), not detection.
  .generate_seir_path <- function(epi_time, X, id, simParam) {
    Tinf   <- epi_time
    Tsym   <- Tinf + rgamma(1L, simParam$LP_shape,
                            scale = simParam$LP_scale * X$lat[[id]])
    Tdeath <- Tsym + rgamma(1L, simParam$RP_shape,
                            scale = simParam$RP_scale * X$tol[[id]])
    list("E", Tinf, Tsym, Tdeath)
  }
  
  # The epidemic is still active while any individual is exposed (E, will become
  # infectious) or infectious (I).
  .get_seir_infectives <- function(X) {
    X[, sum(status %in% c("E", "I"))]
  }
  
  # Next non-infection event. The pending transition depends on the current
  # status: E -> I fires at Tsym, I -> R fires at Tdeath. (epi_time kept for
  # signature parity; status already disambiguates.)
  .next_seir_ni_event <- function(X, epi_time) {
    as.list(X[, .(.I,
                  Tmin = fcase(status == "E", Tsym,
                               status == "I", Tdeath,
                               default = Inf))]
            [, list(t_next_event = min(Tmin, na.rm = TRUE),
                    id_next_event = which.min(Tmin))])
  }
  
  X <- .epi_working_table(epop)
  .epi_init_donors(X, .generate_seir_path, simParam)
  
  #TODO: GE was Removed. Add back
  Xgroups <- X |> split(by = "group")
  
  # Each group is independent, so it seems to be about ~33% faster to split
  # them by group, run them separately, then recombine them.
  Y <- map(Xgroups, \(X) {
    # Start epidemic simulation loop ----
    epi_time <- 0.0
    while (.get_seir_infectives(X) > 0L) {
      if (DEBUG) message("time = ", signif(epi_time, 5))
      
      # Calculate infection rate in the group (only I is infectious) ----
      X[, group_inf := simParam$r_beta *
          mean(fifelse(status == "I", inf, 0.0))]
      
      # if S, infection at rate beta SI
      X[, event_rate := fifelse(status == "S", sus * group_inf, 0.0)]
      
      # id and time of next non-infection event
      ni_event      <- .next_seir_ni_event(X, epi_time)
      t_next_event  <- ni_event$t_next_event
      id_next_event <- ni_event$id_next_event
      
      if (is.na(t_next_event)) t_next_event <- Inf
      
      if (DEBUG) message("Next NI event id = ", id_next_event, " at t = ",
                         signif(t_next_event, 5))
      
      # generate random timestep ----
      total_event_rate <- sum(X$event_rate)
      
      dt <- if (total_event_rate > 0.0) {
        rexp(1L, rate = total_event_rate)
      } else {
        Inf
      }
      
      if (DEBUG) message("Total infections event rate = ",
                         signif(total_event_rate, 5))
      
      # check if next event is infection or non-infection ----
      if (epi_time + dt < t_next_event) {
        if (DEBUG) message("next event is infection at t = ", signif(epi_time, 5))
        
        epi_time <- epi_time + dt
        
        # randomly select individual to be infected
        id_next_event <- sample(nrow(X), size = 1L, prob = X$event_rate)
        
        # select the infector, weighted by infectivity (only I is infectious)
        infectives <- X[status == "I", .(.I, status, inf)]
        infd_by <- .safe_sample(x = infectives$I, size = 1L,
                                prob = infectives$inf)
        next_gen <- X$generation[[infd_by]] + 1L
        
        set(X, id_next_event, c("status", "Tinf", "Tsym", "Tdeath"),
            .generate_seir_path(epi_time, X, id_next_event, simParam))
        set(X, id_next_event, c("generation", "infected_by"),
            list(next_gen, infd_by))
        
        if (DEBUG) message("ID ", id_next_event,
                           ": S -> E, infected by ID ", infd_by)
      } else {
        if (DEBUG) message("next event is non-infection at t = ",
                           signif(t_next_event, 5))
        
        epi_time <- t_next_event
        
        status <- X$status[[id_next_event]]
        
        if (status == "E") {
          if (DEBUG) message("ID ", id_next_event, ": E -> I")
          set(X, id_next_event, "status", "I")
        } else if (status == "I") {
          if (DEBUG) message("ID ", id_next_event, ": I -> R")
          set(X, id_next_event, "status", "R")
        } else {
          message("status = ", status)
          print(X[, .(group, donor, status, Tinf, Tsym, Tdeath,
                      group_inf, event_rate)])
          print(X[id_next_event, .(group, donor, status, Tinf, Tsym, Tdeath,
                                   group_inf, event_rate)])
          stop("selected ID ", id_next_event, "... unexpected event!")
          break
        }
      }
      
      if (DEBUG == 2) {
        print(X[, .(group, donor, status, Tinf, Tsym, Tdeath,
                    group_inf, event_rate)])
      }
    }
    X
  }) |> rbindlist()
  
  .finalise_epi_dynamics(epop, Y)
}

#' @title modelSEIDR
#' @description
#' Compartmental model SEIDR - Susceptibility, Exposed (latent), Infectivity
#'  (undetected), Detectability, Recoverability (or Removed). Exposed (E)
#'  individuals are infected but NOT yet infectious; both undetected (I) and
#'  detected (D) individuals are infectious.
#'
#' @param epop A population with an epidemy, an object of class
#' \code{\link{PopEpidemic-class}}
#' @param simParam a \code{\link{SimParamEpidemic}} object
#'
#' @return Returns the same input object of type \code{\link{PopEpidemic-class}}
#'  with the \code{dynamic} field populated with evolution of the epidemy
modelSEIDR <- function(epop, simParam = NULL) {
  DEBUG <- F
  if (is.null(simParam)) simParam = get("SP", envir = .GlobalEnv)
  
  stopifnot("simParam must be a SimParamEpidemic object" =
              is(simParam, "SimParamEpidemic"),
            "epop must be a PopEpidemic object" =
              is(epop, "PopEpidemic"),
            "Not an SEIDR model" =
              toupper(simParam$model) == "SEIDR")
  
  # A susceptible individual's future trajectory is fixed at infection:
  # E (Tinf) -> I (Tinc) -> D (Tsym) -> R (Tdeath).
  #  - latency  E->I uses 'lat' (LP_*)
  #  - detection I->D uses 'det' (DP_*)
  #  - removal   D->R uses 'tol' (RP_*)
  .generate_seidr_path <- function(epi_time, X, id, simParam) {
    Tinf   <- epi_time
    Tinc   <- Tinf + rgamma(1L, simParam$LP_shape,
                            scale = simParam$LP_scale * X$lat[[id]])
    Tsym   <- Tinc + rgamma(1L, simParam$DP_shape,
                            scale = simParam$DP_scale * X$det[[id]])
    Tdeath <- Tsym + rgamma(1L, simParam$RP_shape,
                            scale = simParam$RP_scale * X$tol[[id]])
    list("E", Tinf, Tinc, Tsym, Tdeath)
  }
  
  # Active while anyone is exposed (E, will become infectious) or infectious
  # (I or D).
  .get_seidr_infectives <- function(X) {
    X[, sum(status %in% c("E", "I", "D"))]
  }
  
  # Next non-infection event. Three possible transitions, chosen by status:
  # E -> I at Tinc, I -> D at Tsym, D -> R at Tdeath. (epi_time kept for
  # signature parity; status already disambiguates.)
  .next_seidr_ni_event <- function(X, epi_time) {
    as.list(X[, .(.I,
                  Tmin = fcase(status == "E", Tinc,
                               status == "I", Tsym,
                               status == "D", Tdeath,
                               default = Inf))]
            [, list(t_next_event = min(Tmin, na.rm = TRUE),
                    id_next_event = which.min(Tmin))])
  }
  
  X <- .epi_working_table(epop)
  .epi_init_donors(X, .generate_seidr_path, simParam)
  
  #TODO: GE was Removed. Add back
  Xgroups <- X |> split(by = "group")
  
  # Each group is independent, so it seems to be about ~33% faster to split
  # them by group, run them separately, then recombine them.
  Y <- map(Xgroups, \(X) {
    # Start epidemic simulation loop ----
    epi_time <- 0.0
    while (.get_seidr_infectives(X) > 0L) {
      if (DEBUG) message("time = ", signif(epi_time, 5))
      
      # Calculate infection rate in the group (both I and D are infectious) ----
      X[, group_inf := simParam$r_beta *
          mean(fifelse(status %in% c("I", "D"), inf, 0.0))]
      
      # if S, infection at rate beta SI
      X[, event_rate := fifelse(status == "S", sus * group_inf, 0.0)]
      
      # id and time of next non-infection event
      ni_event      <- .next_seidr_ni_event(X, epi_time)
      t_next_event  <- ni_event$t_next_event
      id_next_event <- ni_event$id_next_event
      
      if (is.na(t_next_event)) t_next_event <- Inf
      
      if (DEBUG) message("Next NI event id = ", id_next_event, " at t = ",
                         signif(t_next_event, 5))
      
      # generate random timestep ----
      total_event_rate <- sum(X$event_rate)
      
      dt <- if (total_event_rate > 0.0) {
        rexp(1L, rate = total_event_rate)
      } else {
        Inf
      }
      
      if (DEBUG) message("Total infections event rate = ",
                         signif(total_event_rate, 5))
      
      # check if next event is infection or non-infection ----
      if (epi_time + dt < t_next_event) {
        if (DEBUG) message("next event is infection at t = ", signif(epi_time, 5))
        
        epi_time <- epi_time + dt
        
        # randomly select individual to be infected
        id_next_event <- sample(nrow(X), size = 1L, prob = X$event_rate)
        
        # select the infector, weighted by infectivity (from I or D)
        infectives <- X[status %in% c("I", "D"), .(.I, status, inf)]
        infd_by <- .safe_sample(x = infectives$I, size = 1L,
                                prob = infectives$inf)
        next_gen <- X$generation[[infd_by]] + 1L
        
        set(X, id_next_event, c("status", "Tinf", "Tinc", "Tsym", "Tdeath"),
            .generate_seidr_path(epi_time, X, id_next_event, simParam))
        set(X, id_next_event, c("generation", "infected_by"),
            list(next_gen, infd_by))
        
        if (DEBUG) message("ID ", id_next_event,
                           ": S -> E, infected by ID ", infd_by)
      } else {
        if (DEBUG) message("next event is non-infection at t = ",
                           signif(t_next_event, 5))
        
        epi_time <- t_next_event
        
        status <- X$status[[id_next_event]]
        
        if (status == "E") {
          if (DEBUG) message("ID ", id_next_event, ": E -> I")
          set(X, id_next_event, "status", "I")
        } else if (status == "I") {
          if (DEBUG) message("ID ", id_next_event, ": I -> D")
          set(X, id_next_event, "status", "D")
        } else if (status == "D") {
          if (DEBUG) message("ID ", id_next_event, ": D -> R")
          set(X, id_next_event, "status", "R")
        } else {
          message("status = ", status)
          print(X[, .(group, donor, status, Tinf, Tinc, Tsym, Tdeath,
                      group_inf, event_rate)])
          print(X[id_next_event, .(group, donor, status, Tinf, Tinc, Tsym,
                                   Tdeath, group_inf, event_rate)])
          stop("selected ID ", id_next_event, "... unexpected event!")
          break
        }
      }
      
      if (DEBUG == 2) {
        print(X[, .(group, donor, status, Tinf, Tinc, Tsym, Tdeath,
                    group_inf, event_rate)])
      }
    }
    X
  }) |> rbindlist()
  
  .finalise_epi_dynamics(epop, Y)
}
