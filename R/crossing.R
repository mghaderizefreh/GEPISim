# makeCross override -------------------------------------------------------

#' @title Make designed crosses
#'
#' @description
#' Makes crosses within a population using a user supplied
#' crossing plan.
#'
#' @param pop an object of \code{\link{EPop-class}}
#' @param crossPlan a matrix with two columns representing female and male 
#' parents. Either integers for the position in population or character strings
#' for the IDs.
#' @param nProgeny number of progeny per cross
#' @param simParam an object of \code{\link{SimParamEpidemy}}
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @return Returns an object of \code{\link{EPop-class}}
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParamEpidemy$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newEpidemy(founderPop, simParam=SP)
#'
#' #Cross individual 1 with individual 10
#' crossPlan = matrix(c(1,10), nrow=1, ncol=2)
#' pop2 = makeCross(pop, crossPlan, simParam=SP)
#'
#' @export
setMethod(
  "makeCross",
  signature(pop = "EPop"),
  function(pop, crossPlan, nProgeny=1, simParam=NULL, nThreads=NULL){
    # 1. Run the standard makeCross
    pop <- callNextMethod(pop, crossPlan, nProgeny=nProgeny, simParam=simParam, 
                          nThreads=nThreads)
    asEPop(pop,simParam = simParam)
  }
)


