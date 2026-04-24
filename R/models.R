#' @title run
#' @description
#' General wrapper for all compartmental model
#' @param epop A population with an epidemy, an object of class
#' \code{\link{EPop-class}}
#' @param simParam a \code{\link{SimParamEpidemy}} object
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
#' \code{\link{EPop-class}}
#' @param simParam a \code{\link{SimParamEpidemy}} object
#'
#' @return Returns the same input object of type \code{\link{EPop-class}} with
#'  the \code{dynamic} field populated with evolution of the epidemy
modelSIR <- function(epop, simParam = NULL) {
  DEBUG <- F
  if(is.null(simParam)) simParam = get("SP",envir=.GlobalEnv)

  stopifnot("simParam must be a SimParamEpidemy object"=
              is(simParam, "SimParamEpidemy"),
            "epop must be a EPop object"=
              is(epop, "EPop"),
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

  .safe_sample <- function(x, ...) {
    if (length(x) > 1L) {
      sample(x, ...)
    } else if (length(x) <= 1L) {
      x
    }
  }
  # Find the time of the next non-infection event, and the id of the individual
  .next_sir_ni_event <- function(X, epi_time) {
    as.list(X[, .(.I,
                  Tmin = fifelse(Tdeath > epi_time, Tdeath, Inf))]
            [, list(t_next_event = min(Tmin, na.rm = TRUE),
                    id_next_event = which.min(Tmin))])
  }

  #TODO: Quickest way to reuse Jamie's code but maybe not memory friendly
  X <- epop@pheno |> as.data.table() |> cbind(epop@dynamics)

  purrr::walk(X[, .I[donor == 1L]], \(i) {
    data.table::set(X, i, c("status", simParam$timings),
        .generate_sir_path(0.0, X, i, simParam))
  })

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

  final_t <- max(Y$Tdeath, na.rm = TRUE) |> signif(5)

  message(sprintf("- Final t = %f, values are:",final_t),
          paste(capture.output(table(Y$status)), collapse = ", "), "\n")

  # tidy up X
  Y[, c("group_inf", "event_rate") := NULL]
  Y[, parasites := !is.na(Tinf)] # what is this column

  # fix generation
  # X[status == "R" & donor == 0L, generation := 2L]

  epop@dynamics <- Y# rbind(epop@dynamics, Y, fill = TRUE)

  epop
}