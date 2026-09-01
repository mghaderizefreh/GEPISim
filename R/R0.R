#' Get R0
#'
#' Calculate R0 from the output of a simulated epidemic
#'
#' @param ePop A \code{\link{PopEpidemic-class}} object with infection
#'  generations.
#' @param simParam simulation parameter of type \code{\link{SimParamEpidemic}}
#' @param method The method to use. 1 for the ratio of Gen 1 to Gen 2
#' infectives or 2 based on García-Ballesteros et al 202x
#'
#' @returns An estimated R0
#' @export

get_R0 <- function(ePop, simParam = NULL, method = 1) {
  stopifnot("Unknown method passed. Must be 1 or 2" = method%in%c(1,2))

  if (method == 1){
    x <- ePop@dynamics[, table(purrr::map_int(generation, data.table::first))]

    R0 <- if (length(x) > 1) x[[2]] / x[[1]] else 0

  } else if (method == 2){
    if(is.null(simParam)){
      simParam = get("SP",envir=.GlobalEnv)
    }
    if (simParam$model == "SIR"){# for now!
      R0 <-  (simParam$r_beta/simParam$removal_period) *
        mean(ePop@pheno[ePop@dynamics$indCases == 0,'sus']) *
        mean(ePop@pheno[ePop@dynamics$indCases == 1, 'inf'] /
               ePop@pheno[ePop@dynamics$indCases == 1, 'tol'])
    }

  } else{
    R0 <- NULL
  }
  return (R0)

}

