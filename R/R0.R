#' Get R0
#'
#' Calculate R0 given the model and phenotype values
#'
#' @param ePop A \code{\link{PopEpidemic-class}} object with infection
#'  generations.
#' @param simParam simulation parameter of type \code{\link{SimParamEpidemic}}
#'
#' @returns R0 as the basic reproductive ratio
#' @export
get_R0 <- function(ePop, simParam = NULL, method = 1) {
  stopifnot("Unknown method passed. Must be 1 or 2" = method%in%c(1,2))

  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }

  # here we ignore the groups
  if (simParam$model %in% c("SIR", "SEIR")){#
    R0s <- vapply(
      1:ePop@nInd,
      \(i) (
        simParam$r_beta *
          mean(ePop@pheno[-i, 'sus']) *
          ePop@pheno[i, 'inf'] *
          simParam$removal_period * ePop@pheno[i, 'tol']
      ), numeric(1))
    return(median(R0))
  } else if (simParam$model %in%  c("SIDR", "SEIDR")){
    R0s <- vapply(
      1:ePop@nInd,
      \(i) (
        simParm@r_beta *
          mean(ePop@pheno[-i, 'sus']) *
          ePop@pheno[i, 'inf'] *
          (
            simParam$removal_period * ePop@pheno[i, 'tol'] +
              simParam$detection_period * ePop@pheno[i, 'det']
          )
      ), numeric(1))
    return(median(R0))
  } else{
    # not-handled model
    warning(paste("Modle", simParam$model,"is not allowed"))
    return(NA)
  }
}

