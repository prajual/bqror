#' @name  data25j3
#' @title Data containing 500 observations generated from the quantile
#' ordinal model with 3 outcomes and \eqn{p = 0.25} (i.e., 25th quantile)
#'
#' @details
#'
#' This data contains 500 observations generated from the quantile
#' ordinal model with 3 outcomes at the 25th quantile (i.e., \eqn{p = 0.25}).
#' The model specifications for generating the data are as follows: \eqn{\beta = (-4, 6, 5)}, $X ~ Unif(0, 1)$, and
#' \eqn{\epsilon} ~ AL(\eqn{0, \sigma = 1, p = 0.25}).
#'
#' The errors are generated from the asymmetric Laplace distribution
#' using the normal exponential mixture formulation. The cut-points \eqn{(0, 3)}
#' are used to categorize the continuous values into 3 ordinal outcomes.
#'
#' @docType data
#'
#' @usage data(data25j3)
#'
#' @return Returns a list with components
#' \itemize{
#' \item{\code{x}: }{a matrix of covariates.}
#' \item{\code{y}: }{a matrix of ordinal outcomes.}
#' }
#'
#' @references Kozumi, H., and Kobayashi, G. (2011). “Gibbs Sampling Methods for Bayesian Quantile Regression.”
#'  Journal of Statistical Computation and Simulation, 81(11), 1565–1578. DOI: 10.1080/00949655.2010.496117
#'
#' Yu, K., and Zhang, J. (2005). “A Three-Parameter Asymmetric
#' Laplace Distribution.” Communications in Statistics - Theory and Methods, 34(9-10), 1867-1879. DOI: 10.1080/03610920500199018
#'
#' @keywords datasets
#'
#' @seealso \link[MASS]{mvrnorm}, Asymmetric Laplace Distribution
#'
NULL

#' @name  data50j3
#' @title Data containing 500 observations generated from the quantile
#' ordinal model with 3 outcomes and \eqn{p = 0.5} (i.e., 50th quantile)
#'
#' @details
#'
#' This data contains 500 observations generated from the quantile
#' ordinal model with 3 outcomes at the 50th quantile (i.e., \eqn{p = 0.5}).
#' The model specifications for generating the data are as follows: \eqn{\beta = (-4, 6, 5)}, $X ~ Unif(0, 1)$, and
#' \eqn{\epsilon} ~ AL(\eqn{0, \sigma = 1, p = 0.5}).
#'
#' The errors are generated from the asymmetric Laplace distribution
#' using the normal exponential mixture formulation. The cut-points \eqn{(0, 3)}
#' are used to categorize the continuous values into 3 ordinal outcomes.
#'
#' @docType data
#' @usage data(data50j3)
#'
#' @return Returns a list with components
#' \itemize{
#' \item{\code{x}: }{a matrix of covariates.}
#' \item{\code{y}: }{a matrix of ordinal outcomes.}
#' }
#'
#' @references Kozumi, H., and Kobayashi, G. (2011). “Gibbs Sampling Methods for Bayesian Quantile Regression.”
#'  Journal of Statistical Computation and Simulation, 81(11), 1565–1578. DOI: 10.1080/00949655.2010.496117
#'
#' Yu, K., and Zhang, J. (2005). “A Three-Parameter Asymmetric
#' Laplace Distribution.” Communications in Statistics - Theory and Methods, 34(9-10), 1867-1879. DOI: 10.1080/03610920500199018
#'
#' @keywords datasets
#'
#' @seealso \link[MASS]{mvrnorm}, Asymmetric Laplace Distribution
#'
NULL

#' @name  data75j3
#' @title Data containing 500 observations generated from the quantile
#' ordinal model with 3 outcomes and \eqn{p = 0.75} (i.e., 75th quantile)
#'
#' @details
#'
#' This data contains 500 observations generated from the quantile
#' ordinal model with 3 outcomes at the 75th quantile (i.e., \eqn{p = 0.75}).
#' The model specifications for generating the data are as follows: \eqn{\beta = (-4, 6, 5)}, $X ~ Unif(0, 1)$, and
#' \eqn{\epsilon} ~ AL(\eqn{0, \sigma = 1, p = 0.75}).
#'
#' The errors are generated from the asymmetric Laplace distribution
#' using the normal exponential mixture formulation. The cut-points \eqn{(0, 3)}
#' are used to categorize the continuous values into 3 ordinal outcomes.
#'
#' @docType data
#' @usage data(data75j3)
#'
#' @return Returns a list with components
#' \itemize{
#' \item{\code{x}: }{a matrix of covariates.}
#' \item{\code{y}: }{a matrix of ordinal outcomes.}
#' }
#'
#' @references Kozumi, H., and Kobayashi, G. (2011). “Gibbs Sampling Methods for Bayesian Quantile Regression.”
#'  Journal of Statistical Computation and Simulation, 81(11), 1565–1578. DOI: 10.1080/00949655.2010.496117
#'
#' Yu, K., and Zhang, J. (2005). “A Three-Parameter Asymmetric
#' Laplace Distribution.” Communications in Statistics - Theory and Methods, 34(9-10), 1867-1879. DOI: 10.1080/03610920500199018
#'
#' @keywords datasets
#'
#' @seealso \link[MASS]{mvrnorm}, Asymmetric Laplace Distribution
#'
NULL

#' @name  data25j4
#' @title Data containing 500 observations generated from the quantile
#' ordinal model with 4 outcomes and \eqn{p = 0.25} (i.e., 25th quantile)
#'
#' @details
#'
#' This data contains 500 observations generated from the quantile
#' ordinal model with more than 3 outcomes at the 25th quantile (i.e., \eqn{p = 0.25}).
#' The model specifications for generating the data are as follows: \eqn{\beta = (-4, 5, 6)}, $X ~ Unif(0, 1)$, and
#' \eqn{\epsilon} ~ AL(\eqn{0, \sigma = 1, p = 0.25}).
#'
#' The errors are generated from the asymmetric Laplace distribution
#' using the normal exponential mixture formulation. The cut-points \eqn{(0, 2, 4)}
#' are used to categorize the continuous values into 4 ordinal outcomes.
#'
#' @docType data
#' @usage data(data25j4)
#'
#' @return Returns a list with components
#' \itemize{
#' \item{\code{x}: }{a matrix of covariates.}
#' \item{\code{y}: }{a matrix of ordinal outcomes.}
#' }
#'
#' @references Kozumi, H., and Kobayashi, G. (2011). “Gibbs Sampling Methods for Bayesian Quantile Regression.”
#'  Journal of Statistical Computation and Simulation, 81(11), 1565–1578. DOI: 10.1080/00949655.2010.496117
#'
#' Yu, K., and Zhang, J. (2005). “A Three-Parameter Asymmetric
#' Laplace Distribution.” Communications in Statistics - Theory and Methods, 34(9-10), 1867-1879. DOI: 10.1080/03610920500199018
#'
#' @keywords datasets
#'
#' @seealso \link[MASS]{mvrnorm}, Asymmetric Laplace Distribution
#'
NULL

#' @name  data50j4
#' @title Data containing 500 observations generated from the quantile
#' ordinal model with 4 outcomes and \eqn{p = 0.5} (i.e., 50th quantile)
#'
#' @details
#'
#' This data contains 500 observations generated from the quantile
#' ordinal model with more than 3 outcomes at the 50th quantile (i.e., \eqn{p = 0.5}).
#' The model specifications for generating the data are as follows: \eqn{\beta = (-4, 5, 6)}, $X ~ Unif(0, 1)$, and
#' \eqn{\epsilon} ~ AL(\eqn{0, \sigma = 1, p = 0.5}).
#'
#' The errors are generated from the asymmetric Laplace distribution
#' using the normal exponential mixture formulation. The cut-points \eqn{(0, 2, 4)}
#' are used to categorize the continuous values into 4 ordinal outcomes.
#'
#' @docType data
#' @usage data(data50j4)
#'
#' @return Returns a list with components
#' \itemize{
#' \item{\code{x}: }{a matrix of covariates.}
#' \item{\code{y}: }{a matrix of ordinal outcomes.}
#' }
#'
#' @references Kozumi, H., and Kobayashi, G. (2011). “Gibbs Sampling Methods for Bayesian Quantile Regression.”
#'  Journal of Statistical Computation and Simulation, 81(11), 1565–1578. DOI: 10.1080/00949655.2010.496117
#'
#'  Yu, K., and Zhang, J. (2005). “A Three-Parameter Asymmetric
#'  Laplace Distribution.” Communications in Statistics - Theory and Methods, 34(9-10), 1867-1879. DOI: 10.1080/03610920500199018
#'
#' @keywords datasets
#'
#' @seealso \link[MASS]{mvrnorm}, Asymmetric Laplace Distribution
#'
NULL

#' @name  data75j4
#' @title Data containing 500 observations generated from the quantile
#' ordinal model with 4 outcomes and \eqn{p = 0.75} (i.e., 75th quantile)
#'
#' @details
#'
#' This data contains 500 observations generated from the quantile
#' ordinal model with more than 3 outcomes at the 75th quantile (i.e., \eqn{p = 0.75}).
#' The model specifications for generating the data are as follows: \eqn{\beta = (-4, 5, 6)}, $X ~ Unif(0, 1)$, and
#' \eqn{\epsilon} ~ AL(\eqn{0, \sigma = 1, p = 0.75}).
#'
#' The errors are generated from the asymmetric Laplace distribution
#' using the normal exponential mixture formulation. The cut-points \eqn{(0, 2, 4)}
#' are used to categorize the continuous values into 4 ordinal outcomes.
#'
#' @docType data
#' @usage data(data75j4)
#'
#' @return Returns a list with components
#' \itemize{
#' \item{\code{x}: }{a matrix of covariates.}
#' \item{\code{y}: }{a matrix of ordinal outcomes.}
#' }
#'
#' @references Kozumi, H., and Kobayashi, G. (2011). “Gibbs Sampling Methods for Bayesian Quantile Regression.”
#'  Journal of Statistical Computation and Simulation, 81(11), 1565–1578. DOI: 10.1080/00949655.2010.496117
#'
#' Yu, K., and Zhang, J. (2005). “A Three-Parameter Asymmetric
#' Laplace Distribution.” Communications in Statistics - Theory and Methods, 34(9-10), 1867-1879. DOI: 10.1080/03610920500199018
#'
#' @keywords datasets
#'
#' @seealso \link[MASS]{mvrnorm}, Asymmetric Laplace Distribution
#'
NULL

#' @name  Educational_Attainment
#' @title The data consists of information on educational attainment and
#' other variables for 3,923 individuals and is taken from the National
#' Longitudinal Study of Youth (NLSY, 1979) survey.
#'
#' @details
#'
#' The data is taken from the National Longitudinal Study of Youth (NLSY, 1979)
#' survey and corresponds to 3,923 individuals. The objective is to study the
#' effect of family background, individual and school level variable on the
#' quantiles of educational attainment. The dependent variable i.e. the educational
#' degree, has four categories given as less than high school, high school degree,
#' some college or associate's degree and college or graduate degree. The independent
#' variables include intercept, square root of family income, mother's education,
#' father's education, mother's working status, gender, race, and whether the youth
#' lived in an urban area at the age o, and indicator variables to control for age-cohort
#' effects.
#'
#' @docType data
#'
#' @usage data(Educational_Attainment)
#'
#' @return Returns data with components
#' \itemize{
#' \item{\code{mother_work}: }{Indicator for working female at age of 14.}
#' \item{\code{urban}: }{Indicator for the youth living in urban are at age of 14.}
#' \item{\code{south}: }{Indicator for the youth living in South at age of 14.}
#' \item{\code{father_educ}: }{Years of individual's father education.}
#' \item{\code{mother_educ}: }{Years of individual's mother education.}
#' \item{\code{fam_income}: }{Family income of the household in $1000.}
#' \item{\code{female}: }{Indicator for individual's gender.}
#' \item{\code{black}: }{Indicator for black race.}
#' \item{\code{age_cohort_2}: }{Indicator for age of 15.}
#' \item{\code{age_cohort_3}: }{Indicator for age of 16.}
#' \item{\code{age_cohort_4}: }{Indicator for age of 17.}
#' \item{\code{dep_edu_level}: }{matrix of ordinal outcomes.}
#' }
#'
#' @references Rahman, M. A. (2016). “Bayesian
#' Quantile Regression for Ordinal Models.”
#' Bayesian Analysis, 11(1): 1-24. DOI: 10.1214/15-BA939
#'
#' Jeliazkov, I., Graves, J., and Kutzbach, M. (2008). “Fitting and Comparison of Models
#' for Multivariate Ordinal Outcomes.” Advances in Econometrics: Bayesian Econometrics,
#' 23: 115–156. DOI: 10.1016/S0731-9053(08)23004-5
#'
#' Jeliazkov, I., and Rahman, M. A. (2012). “Binary and Ordinal Data Analysis
#' in Economics: Modeling and Estimation” in Mathematical Modeling with Multidisciplinary
#' Applications, edited by X.S. Yang, 123-150. John Wiley & Sons Inc, Hoboken, New Jersey. DOI: 10.1002/9781118462706.ch6
#'
#' @keywords datasets
#'
#' @seealso \href{https://www.bls.gov/nls/nlsy97.htm}{Survey Process}.
#'
NULL

#' @name  Policy_Opinion
#' @title The data consists of public opinion on raising federal income taxes on
#' the rich and a host of other variables for 1,164 individuals and is taken from
#' the 2010-2012 American National Election Studies (ANES) on the Evaluation of
#' Government and Society Study I (EGSS 1)
#'
#' @details
#'
#' The data consists of 1,164 observations taken from the 2010-2012 American National Election
#' Studies (ANES) on the Evaluations of Government and Society Study 1 (EGSS 1). The objective
#' is to analyze public opinion on the proposal to raise federal income taxes for couples (individuals)
#' earning more than $250,000 ($200,000) per year. The responses were recorded as oppose, neither
#' favor nor oppose, or favor the tax increase and forms the dependent variable in the study. The
#' independent variables include indicator variables (or dummy) for employment, income above
#' $75,000, bachelor's and post-bachelor's degree, computer ownership, cellphone ownership, and white race.
#'
#' @docType data
#'
#' @usage data(Policy_Opinion)
#'
#' @return Returns data with components
#' \itemize{
#' \item{\code{Intercept}: }{column of ones.}
#' \item{\code{AgeCat}: }{Indicator for age Category.}
#' \item{\code{IncomeCat}: }{Indicator for household income > $75,000.}
#' \item{\code{Bachelors}: }{Individual's highest degree in Bachelors.}
#' \item{\code{Post.Bachelors}: }{Indicator for highest degree in Masters, Professional or Doctorate.}
#' \item{\code{numComputers}: }{Indicator for computer ownership by individual or household.}
#' \item{\code{CellPhone}: }{Indicator for cellphone ownership by individual or household.}
#' \item{\code{White}: }{Indicator for white race.}
#' \item{\code{y}: }{matrix of ordinal outcomes.}
#' }
#'
#' @references Rahman, M. A. (2016). “Bayesian
#' Quantile Regression for Ordinal Models.”
#' Bayesian Analysis, 11(1): 1-24. DOI: 10.1214/15-BA939
#'
#' Jeliazkov, I., Graves, J., and Kutzbach, M. (2008). “Fitting and Comparison of Models
#' for Multivariate Ordinal Outcomes.” Advances in Econometrics: Bayesian Econometrics,
#' 23: 115–156. DOI: 10.1016/S0731-9053(08)23004-5
#'
#' @keywords datasets
#'
#' @seealso \href{https://electionstudies.org/data-center/}{ANES},  \href{https://georgewbush-whitehouse.archives.gov/cea/progrowth.html}{Tax Policy}
#'
NULL

