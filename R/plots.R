#' Plot Model time series
#'
#' @param ePop an \code{\link{PopEpidemic-class}} object
#' @param simParam simulation parameter of type \code{\link{SimParamEpidemic}}
#' @param trim The percentage of x-axis data to be kept. If the whole x-axis is
#'  shown then it will be difficult to see the early dynamics.

#' @returns A plot of the epidemic
#' @export
plot_model <- function(ePop, simParam, trim = 0.95) {
  plt <- switch(simParam$model,
                "SIS" = plot_SIS(ePop, simParam, trim),
                "SIR" = plot_SIR(ePop, simParam, trim),
                "SEIR" = plot_SEIR(ePop, simParam, trim),
                "SIDR" = plot_SIDR(ePop, simParam, trim),
                "SEIDR" = plot_SEIDR(ePop, simParam, trim))

  plt
}

#' Plot Model time series for SEIDR
#'
#' @param ePop an \code{\link{PopEpidemic-class}} object
#' @param simParam simulation parameter of type \code{\link{SimParamEpidemic}}
#' @param trim The percentage of x-axis data to be kept. If the whole x-axis is
#'  shown then it will be difficult to see the early dynamics.
#'
#' @returns A plot of the epidemic
plot_SEIDR <- function(ePop, simParam, trim) {
  make_time_series_seidr <- function(popn, params, tmax = tmax) {
    compartments <- simParam$compartments |> unlist()

    popn2 <- popn[, .(Tinf, Tinc, Tsign, Tdeath)]
    popn2[is.na(Tsign), Tsign := Tdeath]
    popn2[is.na(Tinc),  Tinc  := Tsign]
    popn2[is.na(Tinf),  Tinf  := Tinc]

    N <- popn2[, .N]

    X <- rbind(
      data.table::data.table(
        event = factor("start", c("start", "infection", "incubation",
                                  "detection", "removal", "end")), time = 0.0),
      data.table::data.table(event = "infection",  time = popn2$Tinf),
      data.table::data.table(event = "incubation", time = popn2$Tinc),
      data.table::data.table(event = "detection",  time = popn2$Tsign),
      data.table::data.table(event = "removal",    time = popn2$Tdeath),
      data.table::data.table(event = "end", time = tmax))

    X <- X[is.finite(time)]

    setorder(X, time)
    #                                                S    E    I    D    R
    X[event == "start",      (compartments) := list(+N,  +0L, +0L, +0L, +0L)]
    X[event == "infection",  (compartments) := list(-1L, +1L, +0L, +0L, +0L)]
    X[event == "incubation", (compartments) := list(+0L, -1L, +1L, +0L, +0L)]
    X[event == "detection",  (compartments) := list(+0L, +0L, -1L, +1L, +0L)]
    X[event == "removal",    (compartments) := list(+0L, +0L, +0L, -1L, +1L)]
    X[event == "end",        (compartments) := list(+0L, +0L, +0L, +0L, +0L)]

    X[, event := NULL]
    X[, S := cumsum(S)]
    X[, E := cumsum(E)]
    X[, I := cumsum(I)]
    X[, D := cumsum(D)]
    X[, R := cumsum(R)]
    X[, ID := I + D]

    rbind(X[time == 0][.N],
          X[time > 0 & time <= tmax]) |>
      data.table::melt("time")
  }

  N <- ePop@nInd

  tmax <- quantile(ePop@dynamics$Tdeath, na.rm = TRUE, prob = trim)[[1]]

  events <- make_time_series_seidr(ePop@dynamics, simParam, tmax)

  ggplot2::ggplot(events) +
    ggplot2::aes(x = time, y = value / N, colour = variable) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::scale_colour_manual(
      "Compartments", breaks = c("S", "E", "I", "D", "ID", "R"),
      labels = c("Susceptible (S)", "Exposed (E)", "Undetectable (I)",
                 "Detectable (D)", "Infectious (I+D)", "Removed (R)"),
      values = c("#46F884FF", "#E1DD37FF", "#F05B12FF", "#30123BFF",
                 "#7A0403FF", "#3E9BFEFF")) +
    ggplot2::coord_cartesian(
      xlim = c(0, min(tmax, max(events$time), na.rm = TRUE)), ylim = c(0, 1)) +
    ggplot2::labs(x = "Time (days)", y = "Proportion", title = "SEIDR model") +
    ggplot2::theme_bw()
}

#' Plot Model time series for SIDR
#'
#' @param ePop an \code{\link{PopEpidemic-class}} object
#' @param simParam simulation parameter of type \code{\link{SimParamEpidemic}}
#' @param trim The percentage of x-axis data to be kept. If the whole x-axis is
#'  shown then it will be difficult to see the early dynamics.
#'
#' @returns A plot of the epidemic
plot_SIDR <- function(ePop, simParam, trim) {
  make_time_series_sidr <- function(popn, params, tmax = tmax) {
    compartments <- simParam$compartments |> unlist()

    popn2 <- popn[, .(Tinf, Tsign, Tdeath)]
    popn2[is.na(Tsign), Tsign := Tdeath]
    popn2[is.na(Tinf),  Tinf  := Tsign]

    N <- popn2[, .N]

    X <- rbind(
      data.table::data.table(
        event = factor("start", c("start", "infection", "detection", "removal",
                                  "end")), time = 0.0),
      data.table::data.table(event = "infection", time = popn2$Tinf),
      data.table::data.table(event = "detection", time = popn2$Tsign),
      data.table::data.table(event = "removal",   time = popn2$Tdeath),
      data.table::data.table(event = "end",       time = tmax))

    X <- X[is.finite(time)]

    setorder(X, time)
    #                                               S    I    D    R
    X[event == "start",     (compartments) := list(+N,  +0L, +0L, +0L)]
    X[event == "infection", (compartments) := list(-1L, +1L, +0L, +0L)]
    X[event == "detection", (compartments) := list(+0L, -1L, +1L, +0L)]
    X[event == "removal",   (compartments) := list(+0L, +0L, -1L, +1L)]
    X[event == "end",       (compartments) := list(+0L, +0L, +0L, +0L)]

    X[, event := NULL]
    X[, S := cumsum(S)]
    X[, I := cumsum(I)]
    X[, D := cumsum(D)]
    X[, R := cumsum(R)]
    X[, ID := I + D]

    rbind(X[time == 0][.N],
          X[time > 0 & time <= tmax]) |>
      data.table::melt("time")
  }

  if(is.null(simParam)){
    simParam = get("SP", envir=.GlobalEnv)
  }

  N <- ePop@nInd

  tmax <- quantile(ePop@dynamics$Tdeath, na.rm = TRUE, prob = trim)[[1]]

  events <- make_time_series_sidr(ePop@dynamics, simParam, tmax)

  ggplot2::ggplot(events) +
    ggplot2::aes(x = time, y = value / N, colour = variable) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::scale_colour_manual(
      "Compartments", breaks = c("S", "I", "D", "ID", "R"),
      labels = c("Susceptible (S)", "Undetectable (I)", "Detectable (D)",
                 "Infectious (I+D)", "Removed (R)"),
      values = c("#3E9BFEFF", "#F05B12FF", "#30123BFF", "#7A0403FF",
                 "#46F884FF")) +
    ggplot2::coord_cartesian(
      xlim = c(0, min(tmax, max(events$time), na.rm = TRUE)), ylim = c(0, 1)) +
    ggplot2::labs(x = "Time (days)", y = "Proportion", title = "SIDR model") +
    ggplot2::theme_bw()
}

#' Plot Model time series for SEIR
#'
#' @param ePop an \code{\link{PopEpidemic-class}} object
#' @param simParam simulation parameter of type \code{\link{SimParamEpidemic}}
#'
#' @returns A plot of the epidemic
plot_SEIR <- function(ePop, simParam, trim) {
  make_time_series_seir <- function(popn, params, tmax = tmax) {
    compartments <- simParam$compartments |> unlist()

    popn2 <- popn[, .(Tinf, Tsign, Tdeath)]
    popn2[is.na(Tsign), Tsign := Tdeath]
    popn2[is.na(Tinf),  Tinf  := Tsign]

    N <- popn2[, .N]

    X <- rbind(
      data.table::data.table(
        event = factor("start", c("start", "infection", "incubation", "removal",
                                  "end")), time = 0.0),
      data.table::data.table(event = "infection",  time = popn2$Tinf),
      data.table::data.table(event = "incubation", time = popn2$Tsign),
      data.table::data.table(event = "removal",    time = popn2$Tdeath),
      data.table::data.table(event = "end",        time = tmax))

    X <- X[is.finite(time)]

    setorder(X, time)
    #                                                S    E    I    R
    X[event == "start",      (compartments) := list(+N,  +0L, +0L, +0L)]
    X[event == "infection",  (compartments) := list(-1L, +1L, +0L, +0L)]
    X[event == "incubation", (compartments) := list(+0L, -1L, +1L, +0L)]
    X[event == "removal",    (compartments) := list(+0L, +0L, -1L, +1L)]
    X[event == "end",        (compartments) := list(+0L, +0L, +0L, +0L)]

    X[, event := NULL]
    X[, S := cumsum(S)]
    X[, E := cumsum(E)]
    X[, I := cumsum(I)]
    X[, R := cumsum(R)]

    rbind(X[time == 0][.N],
          X[time > 0 & time <= tmax]) |>
      data.table::melt("time")
  }

  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }

  N <- ePop@nInd

  tmax <- quantile(ePop@dynamics$Tdeath, na.rm = TRUE, prob = trim)[[1]]

  events <- make_time_series_seir(ePop@dynamics, simParam, tmax)

  ggplot2::ggplot(events) +
    ggplot2::aes(x = time, y = value / N, colour = variable) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::scale_colour_manual(
      "Compartments", breaks = c("S", "E", "I", "R"),
      labels = c("Susceptible", "Exposed", "Infectious", "Removed"),
      values = c("#3E9BFEFF", "#E1DD37FF", "#7A0403FF", "#46F884FF")) +
    ggplot2::coord_cartesian(
      xlim = c(0, min(tmax, max(events$time), na.rm = TRUE)), ylim = c(0, 1)) +
    ggplot2::labs(x = "Time (days)", y = "Proportion", title = "SEIR model") +
    ggplot2::theme_bw()
}

#' Plot Model time series for SIR
#'
#' @param ePop an \code{\link{PopEpidemic-class}} object
#' @param simParam simulation parameter of type \code{\link{SimParamEpidemic}}
#' @param trim The percentage of x-axis data to be kept. If the whole x-axis is
#'  shown then it will be difficult to see the early dynamics.
#'
#' @returns A plot of the epidemic
plot_SIR <- function(ePop, simParam, trim) {
  make_time_series_sir <- function(popn, simParam, tmax = tmax) {
    compartments <- simParam$compartments |> unlist()

    popn2 <- popn[, .(Tinf, Tdeath)]
    popn2[is.na(Tinf), Tinf := Tdeath]

    N <- popn2[, .N]

    X <- rbind(
      data.table::data.table(
        event = factor("start", c("start", "infection","removal", "end")),
        time = 0.0),
      data.table::data.table(event = "infection", time = popn2$Tinf),
      data.table::data.table(event = "removal",   time = popn2$Tdeath),
      data.table::data.table(event = "end",       time = tmax))

    X <- X[is.finite(time)]

    setorder(X, time)
    #                                               S    I    R
    X[event == "start",     (compartments) := list(+N,  +0L, +0L)]
    X[event == "infection", (compartments) := list(-1L, +1L, +0L)]
    X[event == "removal",   (compartments) := list(+0L, -1L, +1L)]
    X[event == "end",       (compartments) := list(+0L, +0L, +0L)]

    X[, event := NULL]
    X[, S := cumsum(S)]
    X[, I := cumsum(I)]
    X[, R := cumsum(R)]

    rbind(X[time == 0][.N],
          X[time > 0 & time <= tmax]) |>
      data.table::melt("time")
  }

  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }

  N <- ePop@nInd

  tmax <- quantile(ePop@dynamics$Tdeath, na.rm = TRUE, prob = trim)[[1]]

  events <- make_time_series_sir(ePop@dynamics, simParam, tmax)

  ggplot2::ggplot(events) +
    ggplot2::aes(x = time, y = value / N, colour = variable) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::scale_colour_manual(
      "Compartments", breaks = c("S", "I", "R"),
      labels = c("Susceptible", "Infectious", "Removed"),
      values = c("#3E9BFEFF", "#7A0403FF", "#46F884FF")) +
    ggplot2::coord_cartesian(
      xlim = c(0, min(tmax, max(events$time), na.rm = TRUE)), ylim = c(0, 1))+
    ggplot2::labs(x = "Time (days)", y = "Proportion", title = "SIR model") +
    ggplot2::theme_bw()
}

#' Plot Model time series for SIS
#'
#' @param ePop an \code{\link{PopEpidemic-class}} object
#' @param simParam simulation parameter of type \code{\link{SimParamEpidemic}}
#'
#' @returns A plot of the epidemic
plot_SIS <- function(ePop, simParam) {
  message("Plotting SIS model")

  N <- ePop[sdp == "progeny", .N]
  tmax <- simParam$tmax

  events <- make_time_series_sis(ePop, simParam)

  ggplot(events) +
    aes(x = time,
        y = value / N,
        colour = variable) +
    geom_line(linewidth = 1.2) +
    scale_colour_manual("Compartments",
                        breaks = c("S", "I"),
                        labels = c("Susceptible", "Infectious"),
                        values = c("blue", "red")) +
    coord_cartesian(xlim = c(0, min(tmax, max(events$time), na.rm = TRUE)),
                    ylim = c(0, 1)) +
    labs(x = "Time (days)",
         y = "Proportion",
         title = "SIS model") +
    theme_bw()
}
