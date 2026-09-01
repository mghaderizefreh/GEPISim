#' @rdname SimParamEpidemic
#' @title Epidemic simulation parameter
#'
#' @description container for global epidemic simulation parameter. The users
#' of this object assume it is stored in environment as "SP", if this object
#' is not explicitly passed. This is the practice in AlphaSimR and SimPlyBee
#'
#' @details This documentation shows details specific to \code{SimParamEpidemic}.
#'   We suggest that you also read all the options provided by the AlphaSimR
#'   \code{\link[AlphaSimR]{SimParam}}. Below we show minimal usage cases for
#'   each \code{SimParamEpidemic} function.
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

    # removal_period field ----
    #' @field removal_period the mean time an individual remains in the "I" or
    #'   infectious state before they are removed.
    removal_period = NA_real_,

    # r_beta field ----
    #' @field r_beta mean transmission rate
    r_beta = NA_real_,

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
    #' @param removal_period see \code{\link[GEPISim]{SimParamEpidemic}}
    #'   field \code{removal_period}
    #'
    #' @param r_beta see \code{\link[GEPISim]{SimParamEpidemic}} field
    #'   \code{r_beta}
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
                          removal_period = 10,
                          r_beta = 0.5,
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
      self$removal_period = removal_period
      self$r_beta <- r_beta
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

      invisible(self)
    },

    ### Overriding Traits methods (Public) ----

    #' @description
    #' Randomly assigns eligible QTLs for one or more additive traits.
    #' If simulating more than one trait, all traits will be pleiotropic
    #' with correlated additive effects. This is a wrapper around the function
    #' \code{\link[AlphaSimR]{SimParam}} method \code{addTraitsA} except the
    #' mean is forced to be zero and names are automatically inferred from the
    #' field \code{epi_traits} of the object
    #' \code{\link[GEPISim]{SimParamEpidemic}}
    #'
    #' @param nQtlPerChr number of QTLs per chromosome. Can be a single value or
    #'   nChr values.
    #' @param mean a vector of mean genetic values for the traits. These will be
    #'   ignored as they are forced to be zero.
    #' @param var a vector of desired genetic variances for the traits. If a
    #'   scalar is provided it will be recycled to provide a variance vector of
    #'   length(name)
    #' @param corA a matrix of correlations between additive effects. If NULL is
    #'   passed, it will replaced by the identity matrix
    #' @param gamma should a gamma distribution be used instead of normal. This
    #'   is not recommended for \code{\link[GEPISim]{SimParamEpidemic}}
    #' @param shape the shape parameter for the gamma distribution
    #'   (the rate/scale parameter of the gamma distribution is accounted
    #'   for via the desired level of genetic variance, the var argument). This
    #'   is not recommended for \code{\link[GEPISim]{SimParamEpidemic}}
    #' @param force should the check for a running simulation be
    #' ignored. Only set to TRUE if you know what you are doing.
    #' @param name Name of the epidemiological traits. if NULL is passed, then
    #'   the names will be inferred from the model name, otherwise they should
    #'   match the field \code{epi_traits} of the object
    #'   \code{\link[GEPISim]{SimParamEpidemic}}
    #' @param nThreads number of threads to use if OpenMP is available.
    #' If \code{NULL}, the number is obtained from \code{self$nThreads}.
    #'
    #' @examples
    #' #Create founder haplotypes
    #' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
    #'
    #' #Set simulation parameters
    #' SP1 = SimParamEpidemic$new(founderPop, "SIR")
    #' SP1$addTraitA(name = c('i','r','s'), nThreads = 1L)
    #'
    #'
    addTraitA = function(nQtlPerChr,mean=0,var=1,corA=NULL,
                         gamma=FALSE,shape=1,force=FALSE,name=NULL,
                         nThreads=NULL){
      if (is.null(name)){
        name = unname(self$epi_traits)
      } else{
        stopifnot(
          "Provided names do not match the names from the model. Pass NULL if in doubt"=
                    all(name%in%unname(self$epi_traits)))
      }
      if (length(mean) == 1){
        mean <- rep(mean, length(name))
      }

      if (!all(mean==0)){
        warning("Non zero mean values were passed. They will be forced to zero")
        mean <- rep(0, length(name))
      }

      if (length(var) == 1){
        var <- rep(var, length(name))
      } else {
        stopifnot("var should have length 1 or length(SP$epi_traits)" =
                    length(var) == length(name))
      }

      if (is.null(corA)){
        corA = diag(rep(1,length(name)))
      } else{
        stopifnot(
          "corA should be NULL or square matrix of size length(SP$epi_traits)" =
                    identical(dim(corA), c(length(name), length(name))))
      }

      super$addTraitA(nQtlPerChr, mean = mean, var = var, corA = corA,
                      gamma = gamma, shape = shape, force = force,
                      name = name, nThreads = nThreads)

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
