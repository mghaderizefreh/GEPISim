#' Plot Model time series
#'
#' @param ePop an \code{\link{PopEpidemic-class}} object
#' @param simParam simulation parameter of type \code{\link{SimParamEpidemic}}
#'
#' @returns A plot of the epidemic
#' @export
plot_model <- function(ePop, simParam) {
  plt <- switch(simParam$model,
                "SIS" = plot_SIS(ePop, simParam),
                "SIR" = plot_SIR(ePop, simParam),
                "SEIR" = plot_SEIR(ePop, simParam),
                "SIDR" = plot_SIDR(ePop, simParam),
                "SEIDR" = plot_SEIDR(ePop, simParam))

  plt
}

#' Plot Model time series for SEIDR
#'
#' @param ePop an \code{\link{PopEpidemic-class}} object
#' @param simParam simulation parameter of type \code{\link{SimParamEpidemic}}
#'
#' @returns A plot of the epidemic
plot_SEIDR <- function(ePop, simParam) {
  message("Plotting SEIDR model")

  N <- ePop[sdp == "progeny", .N]
  tmax <- max(simParam$tmax)

  events <- make_time_series_seidr(ePop, simParam)

  ggplot(events) +
    aes(x = time,
        y = value / N,
        colour = variable) +
    geom_line(linewidth = 1.2) +
    scale_colour_manual("Compartments",
                        breaks = c("S", "E", "I", "D", "ID", "R"),
                        labels = c("Susceptible (S)", "Exposed (E)",
                                   "Undetectable (I)", "Detectable (D)",
                                   "Infectious (I+D)", "Removed (R)"),
                        # values = c("blue", "pink", "purple", "darkgreen", "green", "red")) +
                        values = c(viridisLite::viridis(5), "red")) +
    coord_cartesian(xlim = c(0, min(tmax, max(events$time), na.rm = TRUE)),
                    ylim = c(0, 1)) +
    labs(x = "Time (days)",
         y = "Proportion",
         title = "SEIDR model") +
    theme_bw()
}

#' Plot Model time series for SIDR
#'
#' @param ePop an \code{\link{PopEpidemic-class}} object
#' @param simParam simulation parameter of type \code{\link{SimParamEpidemic}}
#'
#' @returns A plot of the epidemic
plot_SIDR <- function(ePop, simParam) {
  message("Plotting SIDR model")

  N <- ePop[sdp == "progeny", .N]
  tmax <- simParam$tmax

  events <- make_time_series_sidr(ePop, simParam)

  ggplot(events) +
    aes(x = time,
        y = value / N,
        colour = variable) +
    geom_line(linewidth = 1.2) +
    scale_colour_manual("Compartments",
                        breaks = c("S", "I", "D", "ID", "R"),
                        labels = c("Susceptible (S)", "Undetectable (I)",
                                   "Detectable (D)", "Infectious (I+D)",
                                   "Removed (R)"),
                        values = c("blue", "mediumpurple", "purple", "red", "green")) +
    coord_cartesian(xlim = c(0, min(tmax, max(events$time), na.rm = TRUE)),
                    ylim = c(0, 1)) +
    labs(x = "Time (days)",
         y = "Proportion",
         title = "SIDR model") +
    theme_bw()
}

#' Plot Model time series for SEIR
#'
#' @param ePop an \code{\link{PopEpidemic-class}} object
#' @param simParam simulation parameter of type \code{\link{SimParamEpidemic}}
#'
#' @returns A plot of the epidemic
plot_SEIR <- function(ePop, simParam) {
  message("Plotting SEIR model")

  N <- ePop[sdp == "progeny", .N]
  tmax <- simParam$tmax

  events <- make_time_series_seir(ePop, simParam)

  ggplot(events) +
    aes(x = time,
        y = value / N,
        colour = variable) +
    geom_line(linewidth = 1.2) +
    scale_colour_manual("Compartments",
                        breaks = c("S", "E", "I", "R"),
                        labels = c("Susceptible", "Exposed", "Infectious", "Removed"),
                        values = c("blue", "pink", "red", "green")) +
    coord_cartesian(xlim = c(0, min(tmax, max(events$time), na.rm = TRUE)),
                    ylim = c(0, 1)) +
    labs(x = "Time (days)",
         y = "Proportion",
         title = "SEIR model") +
    theme_bw()
}

#' Plot Model time series for SIR
#'
#' @param ePop an \code{\link{PopEpidemic-class}} object
#' @param simParam simulation parameter of type \code{\link{SimParamEpidemic}}
#'
#' @returns A plot of the epidemic
plot_SIR <- function(ePop, simParam) {
  make_time_series_sir <- function(popn, simParam) {
    compartments <- c("S", "I", "R")
    #TODO: use Jamie's formula for tmax
    tmax <- max(popn$Tdeath, na.rm = TRUE) / 4

    popn2 <- popn[, .(Tinf, Tdeath)]
    popn2[is.na(Tinf), Tinf := Tdeath]

    N <- popn2[, .N]

    X <- rbind(
      data.table::data.table(event = factor("start", c("start", "infection",
                                                       "removal", "end")),
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

  N <- ePop@dynamics[, .N]
  # I have no idea what this is
  #tmax <- simParam$tmax
  # I will replace it with final_t
  tmax <- max(ePop@dynamics$Tdeath)

  events <- make_time_series_sir(ePop@dynamics, simParam)

  ggplot2::ggplot(events) +
    ggplot2::aes(x = time,
        y = value / N,
        colour = variable) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::scale_colour_manual(
      "Compartments", breaks = c("S", "I", "R"),
      labels = c("Susceptible", "Infectious", "Removed"),
      values = c("blue", "red", "green")) +
    ggplot2::coord_cartesian(
      xlim = c(0, min(tmax, max(events$time), na.rm = TRUE)),
      ylim = c(0, 1)) +
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
