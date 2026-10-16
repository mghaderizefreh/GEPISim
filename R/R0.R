#' Get R0
#'
#' Calculate R0 given the model and phenotype values
#'
#' @param ePop A \code{\link{PopEpidemic-class}} object with infection
#'  generations.
#' @param simParam simulation parameter of type \code{\link{SimParamEpidemic}}
#' @param aggfunc R0 is a population quantity but here we obtain individual R0s.
#'  The aggfunc determines how these individual quantities are aggregated for a
#'  population level R0. As the distribution is usually skewed only median is 
#'  recommended. 
#'
#' @returns R0 as the basic reproductive ratio
#' @export
getR0 <- function(ePop, simParam = NULL, aggfunc = median) {
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }

  # here we ignore the groups
  if (simParam$model %in% c("SIR", "SEIR", "SIS")){#
    R0s <- vapply(
      1:ePop@nInd,
      \(i) (
        simParam$r_beta *
          mean(ePop@pheno[-i, 'sus']) *
          ePop@pheno[i, 'inf'] *
          simParam$removal_period * ePop@pheno[i, 'tol']
      ), numeric(1))
    return(aggfunc(R0s))
  } else if (simParam$model %in%  c("SIDR", "SEIDR")){
    R0s <- vapply(
      1:ePop@nInd,
      \(i) (
        simParam$r_beta *
          mean(ePop@pheno[-i, 'sus']) *
          ePop@pheno[i, 'inf'] *
          (
            simParam$removal_period * ePop@pheno[i, 'tol'] +
              simParam$detection_period * ePop@pheno[i, 'det']
          )
      ), numeric(1))
    return(aggfunc(R0s))
  } else{
    # not-handled model
    warning(paste("Model", simParam$model,"is not allowed"))
    return(NA)
  }
}

