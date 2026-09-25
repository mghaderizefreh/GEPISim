test_that ("plot_model returns a ggplot object for every implemented model", {
  skip_if_not_installed("ggplot2")

  for (model in c("SIR", "SIDR", "SEIR", "SEIDR")) {
    x <- make_test_epidemic(model = model, r_beta = 0)

    simulated <- suppressMessages(runEpidemic(x$epop, simParam = x$sp))
    plot <- plotModel(simulated, simParam = x$sp)

    expect_s3_class(plot, "ggplot")
  }
})

