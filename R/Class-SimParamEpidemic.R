#' @rdname SimParamEpidemic
#' @title Epidemic simulation parameter
#'
#' @description container for global epidemic simulation parameter. The users
#' of this object assume it is stored in environment as "SP", if this object
#' is not explicitly passed.
#'
#' @details This documentation shows details specific to
#' \code{SimParamEpidemic}. It is strongly recommended that you also read the
#' options provided by the AlphaSimR \code{\link[AlphaSimR]{SimParam}}.
#'
#' @export
SimParamEpidemic <- R6::R6Class(
  classname = "SimParamEpidemic",
  inherit = SimParam,

  # Public ----
  public = list(
    # model field ----
    #' @field model the epidemiological model. The model can only be one of the
    #'   following and is case-insensitive
    #'   SIR,
    model = NA_character_,

    # r_beta field ----
    #' @field r_beta mean transmission rate
    r_beta = NA_real_,

    # removal_period field ----
    #' @field removal_period the mean time an individual remains in the "I" or
    #'   infectious state before they are removed.
    removal_period = NA_real_,

    # RP_shape field ----
    #' @field RP_shape
    #' This is the shape of the gamma distribution for survival. A value of 1
    #' gives an exponential distribution. When RP_shape > 1, it indicates that
    #' hazard generally increases with time (removal becomes more likely the
    #' longer infection lasts). This is common for many disease processes where
    #' risk of recovery/death rises after some days.
    RP_shape = NA_real_,

    # RP_scale field ----
    #' @field RP_scale This is not user defined but calculated as
    #'  ` RP_scale = removal_period / RP_shape `
    RP_scale = NA_real_,

    # detection_period field ----
    #' @field detection_period the mean time an individual stays undetected in
    #'   the "I" state before moving to the detected "D" state. Only consumed by
    #'   models with a "D" compartment (e.g. SIDR).
    detection_period = NA_real_,

    # DP_shape field ----
    #' @field DP_shape shape of the gamma distribution for the detection period.
    #'   A value of 1 gives an exponential distribution.
    DP_shape = NA_real_,

    # DP_scale field ----
    #' @field DP_scale not user defined but calculated as
    #'  ` DP_scale = detection_period / DP_shape `
    DP_scale = NA_real_,

    # latent_period field ----
    #' @field latent_period the mean time an individual stays in the exposed
    #'   "E" (latent, non-infectious) state before becoming infectious "I".
    #'   Only consumed by models with an "E" compartment (e.g. SEIR).
    latent_period = NA_real_,

    # LP_shape field ----
    #' @field LP_shape shape of the gamma distribution for the latent period.
    #'   A value of 1 gives an exponential distribution.
    LP_shape = NA_real_,

    # LP_scale field ----
    #' @field LP_scale not user defined but calculated as
    #'  ` LP_scale = latent_period / LP_shape `
    LP_scale = NA_real_,


    # compartment field ----
    #' @field compartments list of compartment (e.g., "S","I","R" for SIR)
    #' This is a list of one-letter characters determined by the model used
    compartments = list(),

    # timings field ----
    #' @field timings list of strings
    #' This is a list of timings determined by the model used
    timings = list(),

    # epi_traits field ----
    #' @field epi_traits list of strings
    #' This is a list of traits to be used in trait definition. Not defined by
    #' the user
    epi_traits = list(),


    #' @description Starts the process of building a new simulation by creating
    #'   a new SimParamEpidemic object and assigning a founder population of
    #'   genomes to this object.
    #'
    #' @param founderPop an object of
    #'   \code{\link[AlphaSimR:Pop-class]{MapPop-class}}
    #'
    #' @param model see \code{\link[GEPISim]{SimParamEpidemic}} field
    #'   \code{model}
    #'
    #' @param r_beta see \code{\link[GEPISim]{SimParamEpidemic}} field
    #'   \code{r_beta}
    #'
    #' @param removal_period see \code{\link[GEPISim]{SimParamEpidemic}}
    #'   field \code{removal_period}
    #'
    #' @param RP_shape \code{\link[GEPISim]{SimParamEpidemic}} field
    #'   \code{RP_shape}
    #'
    #' @param detection_period \code{\link[GEPISim]{SimParamEpidemic}} field
    #'   \code{detection_period}
    #'
    #' @param DP_shape \code{\link[GEPISim]{SimParamEpidemic}} field
    #'   \code{DP_shape}
    #'
    #' @param latent_period \code{\link[GEPISim]{SimParamEpidemic}} field
    #'   \code{latent_period}
    #'
    #' @param LP_shape \code{\link[GEPISim]{SimParamEpidemic}} field
    #'   \code{LP_shape}
    #'
    #' @examples
    #' founderGenomes <- quickHaplo(nInd = 10, nChr = 3, segSites = 10)
    #' SP <- SimParamEpidemic$new(founderGenomes) # default model is "SIR"
    #' \dontshow{SP$nThreads = 1L}
    #'
    #' # there are only a few models to select from
    #' try(SP <- SimParamEpidemic$new(founderGenomes, model = "ABC"))
    initialize = function(founderPop,
                          model = "SIR",
                          r_beta = 0.5,
                          removal_period = 10,
                          RP_shape = 1,
                          detection_period = 10,
                          DP_shape = 1,
                          latent_period = 10,
                          LP_shape = 1){

      model <- toupper(model)
      stopifnot("provided model is not valid" = model %in% private$.validModels)

      super$initialize(founderPop)
      private$.versionGEPISim <- packageDescription("GEPISim")$Version
      self$model <- model
      self$r_beta <- r_beta

      # removal (or recovery) period - from I (or D if it exists) to R
      self$removal_period = removal_period
      self$RP_shape <- RP_shape
      self$RP_scale <- removal_period / RP_shape

      # detection period (used by D-models such as SIDR)
      self$detection_period <- detection_period
      self$DP_shape <- DP_shape
      self$DP_scale <- detection_period / DP_shape

      # latent period (used by E-models such as SEIR)
      self$latent_period <- latent_period
      self$LP_shape <- LP_shape
      self$LP_scale <- latent_period / LP_shape

      self$compartments <- strsplit(model, split = "")[[1]] |> unique() |>
        as.list()
      private$.GetTimingsAndTraitNames()

      # specifying the log-normality
      self$finalizePheno <- function(pheno, pop, simParam = SP, ...){
        pheno <- asLogNormal(pheno)
        return(pheno)
      }

      invisible(self)
    }

  ),
  private = list(
    #### Private ----
    .versionGEPISim = "character",

    .validModels = c("SIR", "SIDR", "SEIR", "SEIDR"),# c("SI","SIS")


    .all_traits = c(s = "sus", i = "inf", l = "lat", d = "det",
                    t = "tol"),
    .GetTimingsAndTraitNames = function(){
      switch(toupper(self$model),
             "SEIDR" = {
               self$epi_traits <- private$.all_traits
               self$timings <- c("Tinf", "Tinc", "Tsign", "Tdeath")
             }, "SIDR" = {
               self$epi_traits <- private$.all_traits[c("s", "i", "d", "t")]
               self$timings <- c("Tinf", "Tsign", "Tdeath")
             }, "SEIR" = {
               self$epi_traits <- private$.all_traits[c("s", "i", "l", "t")]
               self$timings <- c("Tinf", "Tsign", "Tdeath")
             }, "SIR" = {
               self$epi_traits <- private$.all_traits[c("s", "i", "t")]
               self$timings <- c("Tinf", "Tdeath")
             }, "SIS" = {
               self$epi_traits <- private$.all_traits[c("s", "i", "t")]
               self$timings <- c("Tinf", "Tdeath")
             }, "SI" = {
               self$epi_traits <- private$.all_traits[c("s", "i")]
               self$timings <- c("Tinf")
             }, {
               stop("- Unknown model!")
             }
      )
      invisible(self)
    }
  ),

  active = list(
    #' @field version list, versions of AlphaSimR and GEPISim packages used
    #'   to generate this object
    version = function(value){
      if(missing(value)){
        list(
          "AlphaSimR" = private$.version,
          "GEPISim" = private$.versionGEPISim
        )
      } else {
        stop("`$version` is read only", call. = F)
      }
    }
  )
)
