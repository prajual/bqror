#' Bayesian quantile regression for ordinal models
#'
#' @description
#'
#' Package provides functions for estimating Bayesian quantile regression with ordinal outcomes,
#' computing the covariate effects, model comparison measures, and inefficiency factor. The generic
#' ordinal model with 3 or more outcomes (labeled OR1 model) is estimated by a combination of Gibbs
#' sampling and Metropolis-Hastings algorithm. Whereas an ordinal model with exactly 3 outcomes
#' (labeled OR2 model) is estimated using Gibbs sampling only. For each model framework, there is a
#' specific function for estimation. The summary output produces estimates for regression quantiles
#' and two measures of model comparison — log of marginal likelihood and Deviance Information Criterion
#' (DIC). The package also has specific functions for computing the covariate effects and other functions
#' that aids either the estimation or inference in quantile ordinal models
#'
#' @details
#' \deqn{Package: bqror}
#' \deqn{Type: Package}
#' \deqn{Version: 1.4.0}
#' \deqn{License: GPL (>=2)}
#'
#' Package \strong{bqror} provides the following functions:
#'
#' \itemize{
#' \item{For an ordinal model with three or more outcomes:}
#' }
#' \code{\link[bqror]{quantregOR1}}, \code{\link[bqror]{covEffectOR1}},
#' \code{\link[bqror]{logMargLikeOR1}}, \code{\link[bqror]{devianceOR1}},
#' \code{\link[bqror]{qrnegLogLikensumOR1}}, \code{\link[bqror]{infactorOR1}},
#' \code{\link[bqror]{qrminfundtheorem}}, \code{\link[bqror]{drawbetaOR1}},
#' \code{\link[bqror]{drawwOR1}}, \code{\link[bqror]{drawlatentOR1}},
#' \code{\link[bqror]{drawdeltaOR1}}, \code{\link[bqror]{alcdfstd}},
#' \code{\link[bqror]{alcdf}}, \code{\link[bqror]{logLik.bqrorOR1}}
#'
#' \itemize{
#' \item{For an ordinal model with three outcomes:}
#' }
#' \code{\link[bqror]{quantregOR2}}, \code{\link[bqror]{covEffectOR2}},
#' \code{\link[bqror]{logMargLikeOR2}}, \code{\link[bqror]{devianceOR2}},
#' \code{\link[bqror]{qrnegLogLikeOR2}}, \code{\link[bqror]{infactorOR2}},
#' \code{\link[bqror]{drawlatentOR2}}, \code{\link[bqror]{drawbetaOR2}},
#' \code{\link[bqror]{drawsigmaOR2}}, \code{\link[bqror]{drawnuOR2}},
#' \code{\link[bqror]{rndald}}
#'
#' @author
#' Mohammad Arshad Rahman
#'
#' Prajual Maheshwari <prajual1391@gmail.com>
#'
#' @references
#' Rahman, M. A. (2016). “Bayesian
#' Quantile Regression for Ordinal Models.”
#' Bayesian Analysis, 11(1): 1-24. DOI: 10.1214/15-BA939
#'
#' Yu, K., and Moyeed, R. A. (2001). “Bayesian Quantile Regression.” Statistics and
#' Probability Letters, 54(4): 437–447. DOI: 10.1016/S0167-7152(01)00124-9
#'
#' Koenker, R., and Bassett, G. (1978).“Regression Quantiles.” Econometrica,
#' 46(1): 33-50. DOI: 10.2307/1913643
#'
#' Greenberg, E. (2012). “Introduction to Bayesian Econometrics.”
#' Cambridge University Press. Cambridge, DOI: 10.1017/CBO9781139058414
#'
#' @seealso \link[GIGrvg]{rgig}, \link[MASS]{mvrnorm}, \link[MASS]{ginv},
#' \link[truncnorm]{rtruncnorm}, \link[NPflow]{mvnpdf},
#' \link[invgamma]{rinvgamma}, \link[pracma]{mldivide},
#' \link[pracma]{rand}, \link[stats]{qnorm},
#' \link[stats]{rexp}, \link[stats]{rnorm},
#' \link[pracma]{std}, \link[stats]{sd}, \link[stats]{acf},
#' \link[pracma]{Reshape}, \link[progress]{progress_bar},
#' \link[invgamma]{dinvgamma}, \link[stats]{logLik}
#'
#' @docType package
#' @name bqror
NULL


