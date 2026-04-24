# EPop --------------------------------------------------------------------

#' @title Epidemy Population
#'
#' @description
#' A population with an epidemy. Extends \code{\link[AlphaSimR:Pop-class]{Pop}}
#' class from AlphaSimR package to include a new slot and override some of
#' the methods.
#'
#' @param object a 'Pop' object
#' @param x a 'Pop' object
#' @param i index of individuals
#' @param ... additional 'Pop' objects
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
#' open slot available for uses to store extra information about
#' individuals.
#' @slot miscPop a list of any length containing optional meta data for the
#' population (see example in \code{\link[AlphaSimR:newPop]{newPop}}).
#' This list is empty unless information is supplied by the user.
#' Note that the list is emptied every time the population is subsetted or
#' combined because the meta data for old population might not be valid anymore.
#' @name EPop-class
#' @rdname EPop-class
#' @exportClass EPop
setClass("EPop",
         slots=c(dynamics = 'ANY'),
         contains="Pop")
#TODO: setValidity should check for certain characteristics of `dyanimc`

#' @title Create new Epopulation
#'
#' @description
#' Creates an initial \code{\link{EPop-class}} from an object of
#' \code{\link[AlphaSimR:MapPop-class]{MapPop}} or
#' \code{\link[AlphaSimR:NamedMapPop-class]{NamedMapPop}}.
#' This will use AlphaSimR `new` function for
#' \code{\link[AlphaSimR:Pop-class]{Pop}}. Then, it will also embed relevant
#' info from SP in it.
#'
#' @param rawPop an object of \code{\link[AlphaSimR:MapPop-class]{MapPop}} or
#' \code{\link[AlphaSimR:NamedMapPop-class]{NamedMapPop}}
#' @param simParam an object of \code{\link{SimParamEpidemy}}
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from code_simParam$nThreads.
#' @param ... additional arguments used internally (passed to \code{\link[AlphaSimR:Pop-class]{Pop}})
#'
#' @return Returns an object of \code{\link{EPop-class}}
#'
#' @seealso \code{\link[AlphaSimR:Pop-class]{Pop}}, \code{\link{SimParamEpidemy}}
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=2, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParamEpidemy$new(founderPop, model = "SIR")
#' SP$addTraitA(var = c(1,1,1))
#'
#' #Create population
#' pop = newEpidemy(founderPop)
#'
#' @export
newEpidemy = function(rawPop, group_list = NULL, donor_list = NULL,
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
  if (is.null(group_list)){ # if no group provided, assume only 1 group
    group_list <- rep(1L, pop@nInd)
  }
  if (is.null(donor_list)){# if no donor provided, one donor per group
    donor_list <- rep(0, length(group_list))
    groups <- unique(group_list)
    donor_list[sapply(1:length(groups), \(i) (sample(which(
      group_list == groups[i]),1)))] <- 1
  }

  # some checks
            # Validate simParam is SimParamEpidemy
  stopifnot("simParam must be a SimParamEpidemy object"=
              is(simParam, "SimParamEpidemy"),
            #Ensure group_list length matches number of individuals
            "group_list must have equal elements to number of individuals"=
              (length(group_list) == pop@nInd),
            # Ensure donor_list has same length as group_list
            "donor_list must have the same length as group_list"=
              (length(donor_list) == length(group_list)),
            # Ensure donor_list contains only 0s and 1s
            "donor_list must contain only 0 and 1 values"=
              all(donor_list %in% c(0,1))
            )
  # Ensure at least one donor (1) per group
  for(g in unique(group_list)){
    if(sum(donor_list[group_list == g] == 1) < 1){
      stop(sprintf("Each group must have at least one donor (value 1) in donor_list. Group: %s", as.character(g)))
    }
  }
  # add donor_list and group_list to dynamics
  dynamics <- data.table(donor = donor_list, group = group_list)

  # add timing columns
  dynamics[,(simParam$timings) := NA_real_]

  compartments <- strsplit(simParam$model, split = "")[[1]] |> unique()
  # adding 'status' column as factor with levels of the model
  # add group_inf (infection force)
  dynamics[, `:=`(status = factor("S", levels = compartments),
                  group_inf = 0.0)]

  # add columns generation (of infection) and infected_by
  dynamics[donor == 1L, `:=`(generation = 1L, infected_by = 0L)]

  output <- new("EPop",
               nInd=pop@nInd,
               nChr=pop@nChr,
               ploidy=pop@ploidy,
               nLoci=pop@nLoci,
               sex=pop@sex,
               geno=pop@geno,
               id=pop@id,
               iid=pop@iid,
               mother=pop@mother,
               father=pop@father,
               fixEff=rep(1L,pop@nInd),
               misc=list(),
               miscPop=list(),
               nTraits=simParam$nTraits,
               gv=pop@gv,
               gxe=pop@gxe,
               pheno=pop@pheno,
               ebv=matrix(NA_real_,
                          nrow=pop@nInd,
                          ncol=0),
               dynamics = dynamics)
  return(output)
}

#' @describeIn EPop Show epidmey population summary
#' @export
setMethod("show",
          signature(object = "EPop"),
          function (object){
            cat("An object of class",
                classLabel(class(object)), "\n")
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
