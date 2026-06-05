# makeCross override -------------------------------------------------------

#' @title Make designed crosses
#'
#' @description
#' Makes crosses within a population using a user supplied
#' crossing plan.
#'
#' @param pop an object of \code{\link{PopEpidemic-class}}
#' @param crossPlan a matrix with two columns representing female and male
#' parents. Either integers for the position in population or character strings
#' for the IDs.
#' @param nProgeny number of progeny per cross
#' @param simParam an object of \code{\link{SimParamEpidemic}}
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @return Returns an object of \code{\link{PopEpidemic-class}}
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParamEpidemic$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPopEpidemic(founderPop, simParam=SP)
#'
#' #Cross individual 1 with individual 10
#' crossPlan = matrix(c(1,10), nrow=1, ncol=2)
#' pop2 = makeCross(pop, crossPlan, simParam=SP)
#'
#' @export
setMethod(
  "makeCross",
  signature(pop = "PopEpidemic"),
  function(pop, crossPlan, nProgeny=1, simParam=NULL, nThreads=NULL){
    # 1. Run the standard makeCross
    pop <- callNextMethod(pop, crossPlan, nProgeny=nProgeny, simParam=simParam,
                          nThreads=nThreads)
    asPopEpidemic(pop,simParam = simParam)
  }
)

# makeCross2 override -------------------------------------------------------
#' @title Make designed crosses
#'
#' @description
#' Makes crosses between two populations using a user supplied
#' crossing plan.
#'
#' @param females an object of \code{\link{Pop-class}} for female parents.
#' @param males an object of \code{\link{Pop-class}} for male parents.
#' @param crossPlan a matrix with two column representing
#' female and male parents. Either integers for the position in
#' population or character strings for the IDs.
#' @param nProgeny number of progeny per cross
#' @param simParam an object of \code{\link{SimParam}}
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @return Returns an object of \code{\link{Pop-class}}
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#'
#' #Cross individual 1 with individual 10
#' crossPlan = matrix(c(1,10), nrow=1, ncol=2)
#' pop2 = makeCross2(pop, pop, crossPlan, simParam=SP)
#'
#' @export
setMethod(
  "makeCross2",
  signature(females = "PopEpidemic", males = "PopEpidemic"),
  function(females,males,crossPlan,nProgeny=1,simParam=NULL, nThreads=NULL){
    pop <- callNextMethod(females, males, crossPlan, nProgeny=nProgeny,
                          simParam = simParam, nThreads = nThreads)
    asPopEpidemic(pop, simParam = simParam)
  }
)


