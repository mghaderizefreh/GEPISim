# setPheno override -------------------------------------------------------

#' @title Sets phenotypes for epidemiological traits
#'
#' @description
#' Extends \code{\link[AlphaSimR]{setPheno}} for objects of class
#' \code{\link{EPop-class}} to scale multiplicative traits in an epidemiological
#' model such that they have an expected mean of 1 (if desired).
#'
#' @param pop an object of \code{\link{EPop-class}}
#' @param h2 a vector of desired narrow-sense heritabilities for
#' each trait. See details.
#' @param H2 a vector of desired broad-sense heritabilities for
#' each trait. See details.
#' @param varE error (co)variances for traits. See details.
#' @param corE an optional matrix for correlations between errors.
#' See details.
#' @param reps number of replications for phenotype. See details.
#' @param fixEff fixed effect to assign to the population. Used
#' by genomic selection models only.
#' @param p the p-value for the environmental covariate
#' used by GxE traits. If NULL, a value is
#' sampled at random.
#' @param onlyPheno should only the phenotype be returned, see return
#' @param traits an integer vector indicate which traits to set. If NULL,
#' all traits will be set.
#' @param simParam an object of \code{\link{SimParam}}
#' @param stnd how to standardise the exponential traits. Possible values are
#'  'shift' and 'divide'. NULL or other string characters will only exponentiate
#'  the traits without standardising in any way. See details for how each works.
#'
#' @details
#' For description on the parent function see
#' \code{\link[AlphaSimR]{setPheno}}.
#'
#' Shifting transformation:
#'
#' \eqn{P=e^{P-\sigma^2/2}}.
#'
#' Dividing works by transformation
#'
#' \eqn{P=\frac{e^{P}}{\mu(e^{P})}},
#'
#' where \eqn{\mu} is the mean function.
#'
#' \code{NULL} or other character values does the exponentiation trnasformation
#' only
#'
#' \eqn{P=e^{P}}
#' @name setPheno
#' @rdname setPheno
#' @export
setMethod(
  "setPheno",
  signature(pop = "EPop"),
  function(pop, h2 = NULL, H2 = NULL, varE = NULL, corE = NULL,
           reps = 1, fixEff = 1L, p = NULL, onlyPheno = FALSE,
           traits = NULL, simParam = NULL, stnd = NULL) {

    # 1. Run the standard Pop phenotype calculation.
    pop <- callNextMethod(pop, h2, H2, varE, corE, reps, fixEff, p,
                         onlyPheno, traits, simParam)

    # 2. shift by half of the variance and raise to power e
    # TODO: fixed effects are assuming to be not present at the moment

    if(is.null(simParam)) simParam <- get("SP",envir=.GlobalEnv)
    if(is.null(traits)){
      traits <- if (simParam$nTraits > 0L) 1:simParam$nTraits else integer()
    }else{
      traits <- as.integer(traits)
    }
    nTraits <- length(traits)

    varG <- simParam$varG[traits]
    # Calculate varE if using h2 or H2
    if(!is.null(h2)){
      if(length(h2)==1) h2 <- rep(h2, nTraits)
      varA <- simParam$varA[traits]
      varE <- varA[traits]/h2[traits]-varG[traits]
    }else if(!is.null(H2)){
      if(length(H2)==1) H2 <- rep(H2, nTraits)
      varE <- varG[traits]/H2[traits]-varG[traits]
    }else if(!is.null(varE)){
      # varE already given
    }else{
      if(is.matrix(simParam$varE)){
        varE <- simParam$varE[traits, traits]
      }else{
        varE <- simParam$varE[traits]
      }
    }

    # get half of the total phenotypic variance
    varP_half <- (varE + varG) / 2.

    if (is.null(stnd)){
      pop@phen <- exp(pop@pheno)
    }
    if (stnd == 'shift'){
      #shift all the phenotypes and then exponentiate
      shift <- sweep(pop@pheno, 2, varP_half)
      pop@pheno <- exp(shift)
    } else if (stnd == 'divide'){
      # exponentiate then
      popexp <- exp(pop@pheno)
      pop@pheno <- sweep(popexp, 2, colMeans(popexp), `/`)
    }
    else{
      pop@phen <- exp(pop@pheno)
    }
    return(pop)
  })
