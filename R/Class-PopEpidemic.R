# PopEpidemic --------------------------------------------------------------------

#' @title Epidemic Population
#'
#' @description
#' A population with an epidemy. Extends \code{\link[AlphaSimR:Pop-class]{Pop}}
#' class from AlphaSimR package to include a new slot and override some of
#' the methods.
#'
#' @slot id an individual's identifier
#' @slot iid an individual's internal identifier
#' @slot mother the identifier of the individual's mother
#' @slot father the identifier of the individual's father
#' @slot sex sex of individuals: "M" for males, "F" for females,
#' and "H" for hermaphrodites
#' @slot nTraits number of traits
#' @slot gv matrix of genetic values. When using GxE traits,
#' gv reflects gv when p=0.5. Dimensions are nInd by nTraits.
#' @slot pheno matrix of phenotypic values. Dimensions are
#' nInd by nTraits.
#' @slot ebv matrix of estimated breeding values. Dimensions
#' are nInd rows and a variable number of columns.
#' @slot gxe list containing GxE slopes for GxE traits
#' @slot fixEff a fixed effect relating to the phenotype.
#' Used by genomic selection models but otherwise ignored.
#' @slot misc a list whose elements correspond to additional miscellaneous
#' nodes with the items for individuals in the population (see example in
#' \code{\link[AlphaSimR:newPop]{newPop}}) - we support vectors and matrices or
#' objects that have a generic length and subset method.
#' This list is normally empty and exists solely as an
#' open slot available for users to store extra information about
#' individuals.
#' @slot miscPop a list of any length containing optional meta data for the
#' population (see example in \code{\link[AlphaSimR:newPop]{newPop}}).
#' This list is empty unless information is supplied by the user.
#' Note that the list is emptied every time the population is subsetted or
#' combined because the meta data for old population might not be valid anymore.
#' @name PopEpidemic-class
#' @rdname PopEpidemic-class
#' @exportClass PopEpidemic
setClass("PopEpidemic",
         slots=c(dynamics = 'ANY'),
         contains="Pop")

#' @title Create a new population with epidemy (PopEpidemic)
#'
#' @description
#' Creates an initial \code{\link{PopEpidemic-class}} from an object of
#' \code{\link[AlphaSimR:MapPop-class]{MapPop}} or
#' \code{\link[AlphaSimR:NamedMapPop-class]{NamedMapPop}}.
#' This will use AlphaSimR `newPop` function for
#' \code{\link[AlphaSimR:Pop-class]{Pop}}. Then, it will also embed relevant
#' info from SP in it.
#'
#' @param rawPop an object of \code{\link[AlphaSimR:MapPop-class]{MapPop}} or
#' \code{\link[AlphaSimR:NamedMapPop-class]{NamedMapPop}}
#' @param group_list list of groups for the population. Each group is a separate
#'   epidemy. If NULL it is assumed there is only one group
#' @param donor_list list of 0s and 1s indicating susceptible and donor 
#'   individuals, respectively. If any other number is passed this will cause an
#'   error. If NULL, it is assumed there is one donor per group.
#' @param simParam an object of \code{\link{SimParamEpidemic}}
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#' @param ... additional arguments used internally (passed to 
#' \code{\link[AlphaSimR:Pop-class]{newPop}})
#'
#' @return Returns an object of \code{\link{PopEpidemic-class}} with a new slot 
#'  `dynamics`. See \code{\link{initDynamics}} for details
#'
#' @seealso \code{\link[AlphaSimR:Pop-class]{Pop}}, \code{\link{SimParamEpidemic}}
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=2, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParamEpidemic$new(founderPop, model = "SIR")
#' SP$addTraitA(var = c(1,1,1))
#'
#' #Create population
#' pop = newPopEpidemic(founderPop)
#'
#' @export
newPopEpidemic <- function(rawPop, group_list = NULL, donor_list = NULL,
                      simParam=NULL,nThreads=NULL,...){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }
  pop <- AlphaSimR::newPop(rawPop, simParam, nThreads, ...)

  asPopEpidemic(pop, simParam = simParam, group_list = group_list,
         donor_list = donor_list)
}

#' @describeIn PopEpidemic-class Show epidemy population summary
#' @export
setMethod("show",
          signature(object = "PopEpidemic"),
          function (object){
            cat("An object of class",
                class(object), "\n")
            cat("Ploidy:", object@ploidy,"\n")
            cat("Individuals:", object@nInd,"\n")
            cat("Chromosomes:", object@nChr,"\n")
            cat("Loci:", sum(object@nLoci),"\n")
            cat("Traits:", object@nTraits,"\n")
            cat("Model:", paste0(levels(object@dynamics$status),collapse = ''),
                "\n")
            invisible()
          }
)

#' @describeIn PopEpidemic-class initialises the dynamic object
#' @param ePop an \code{\link{PopEpidemic-class}} object
#' @param simParam simulation parameter of type \code{\link{SimParamEpidemic}}
#' @return the same input object ePop with the field `dynamics` initialised.
#' @details
#' The dynamics field is data.table with several columns,
#' - Timing columns: these depend on the model used. They are supposed to be
#'   stored in the simParam and will be extracted from it
#' - status: these also depend on the model used. In general they are the 
#'   letters that define the model, e.g., SIR will have states S, I, and R
#' - generation: Which type of infection is it, i.e., primary, secondary, etc.
#'   Initially set to 1 for all donors
#' - infected_by: the index of the individual that infects a given individual 
#' - group_inf: group infectivity level
#' @export
initDynamics <- function(ePop, simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  # some checks
  stopifnot(
    # Validate simParam is SimParamEpidemic
    "simParam must be a SimParamEpidemic object"=
      is(simParam, "SimParamEpidemic"),
    "object must be an PopEpidemic object"=
      is(ePop, "PopEpidemic")
  )
  
  # add timing columns
  ePop@dynamics[,(simParam$timings) := NA_real_]
  
  compartments <- strsplit(simParam$model, split = "")[[1]] |> unique()
  # adding 'status' column as factor with levels of the model
  # add group_inf (infection force)
  ePop@dynamics[, `:=`(status = factor("S", levels = compartments),
                       group_inf = 0.0)]
  
  # add columns generation (of infection) and infected_by
  ePop@dynamics[donor == 1L, `:=`(generation = 1L, infected_by = 0L)]
  
  ePop
}

#' @describeIn PopEpidemic-class coerces into PopEpidemic-class from Pop-class
#' @param from source object (Pop-class)
#' @param group_list list of groups for the population. Each group is a separate
#'   epidemy. If NULL it is assumed there is only one group
#' @param donor_list list of 0s and 1s indicating susceptible and donor 
#'   individuals, respectively. If any other number is passed this will cause an
#'   error. If NULL, it is assumed there is one donor per group.
#' @param simParam simulation parameter (of type SimParamEpidemic)
#' @export
asPopEpidemic <- function(from, group_list = NULL, donor_list = NULL, simParam = NULL){
  if (is.null(group_list)){ # if no group provided, assume only 1 group
    group_list <- rep(1L, from@nInd)
  }
  if (is.null(donor_list)){# if no donor provided, one donor per group
    donor_list <- rep(0, length(group_list))
    groups <- unique(group_list)
    donor_list[sapply(1:length(groups), \(i) (sample(which(
      group_list == groups[i]),1)))] <- 1
  }
  
  if(is.null(simParam)) simParam = get("SP",envir=.GlobalEnv)
  
  # some checks
  stopifnot(
    # Validate simParam is SimParamEpidemic
    "simParam must be a SimParamEpidemic object"=
      is(simParam, "SimParam"),
    #Ensure group_list length matches number of individuals
    "group_list must have equal elements to number of individuals"=
      (length(group_list) == from@nInd),
    #Ensure donor_list has same length as group_list
    "donor_list must have the same length as group_list"=
     (length(donor_list) == length(group_list)),
    # Ensure donor_list contains only 0s and 1s
    "donor_list must contain only 0 and 1 values"=
      all(donor_list %in% c(0,1))
  )
  # More checks: Ensure at least one donor (1) per group
  for(g in unique(group_list)){
    if(sum(donor_list[group_list == g] == 1) < 1){
      stop(sprintf("Each group must have at least one donor (value 1) in donor_list. Group: %s", as.character(g)))
    }
  }
  
  # add donor_list and group_list to dynamics
  dynamics <- data.table(donor = donor_list, group = group_list)
  
  out <- new(
    "PopEpidemic", nInd=from@nInd, nChr=from@nChr, ploidy=from@ploidy, 
    nLoci=from@nLoci, sex=from@sex, geno=from@geno, id=from@id, iid=from@iid, 
    mother=from@mother, father=from@father, fixEff=rep(1L,from@nInd), 
    misc=list(), miscPop=list(), nTraits=simParam$nTraits, gv=from@gv, 
    gxe=from@gxe, pheno=from@pheno, ebv=matrix(NA_real_,nrow=from@nInd,ncol=0), 
    dynamics = dynamics
  )
  
  initDynamics(out, simParam)
}
