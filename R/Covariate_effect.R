#' Covariate Effect for Bayesian Quantile Regression for Ordinal Models
#'
#' This function estimates the change in probability of different ordinal
#' outcomes due to change in an independent variable, marginalized over
#' the parameters and values of other covariates
#'
#' @param model     outcome of the ODR I (quan_regg3) model.
#' @param y         dependent variable i.e. ordinal outcome values.
#' @param x         covariate matrix of dimension \eqn{(n x k)} including a column of ones with or without column names.
#' @param mod_x     matrix x with suitable modification to an independent variable including a column of ones with or without column names.
#' @param p         quantile level or skewness parameter, p in (0,1).
#'
#' @details
#' Function estimates the covariate effect of the change in an independent
#' variable. The effect is marginalized over other parameters and values of other covariates.
#' This function can only be used for ordinal models with more than 3 outcomes. Function uses
#' the output of ODR I model to compute the covariate effect.
#'
#' @return Returns a list with components:
#' \itemize{
#' \item{\code{avgdiff_prob}: }{a vector with change in predicted
#' probabilities for each outcome category.}
#' }
#'
#' @references Rahman, M. A. (2016). “Bayesian
#' Quantile Regression for Ordinal Models.”
#' Bayesian Analysis, 11(1): 1-24.
#'
#' Jeliazkov, I., Graves, J., and Kutzbach, M. (2008). “Fitting and Comparison of Models
#' for Multivariate Ordinal Outcomes.” Advances in Econometrics: Bayesian Econo-
#' metrics, 23: 115–156.
#'
#' @importFrom "stats" "sd"
#' @examples
#'  set.seed(101)
#'  data("data25j4")
#'  x <- data25j4$x
#'  y <- data25j4$y
#'  p <- 0.25
#'  output <- quan_regg3(y, x, mc = 50, p, 0.1)
#'  mod_x <- x
#'  mod_x[,3] <- mod_x[,3] + 0.02
#'  res <- covariate_effect(output, y, x, mod_x, p)
#'
#'  # avgdiff_prob
#'  #   -4.156752e-03  4.945857e-05 -4.134889e-03  4.085430e-03
#'
#' @export
covariate_effect <- function(model, y, x, mod_x, p) {
    cols <- colnames(x)
    cols1 <- colnames(mod_x)
    names(mod_x) <- NULL
    names(y) <- NULL
    names(x) <- NULL
    x <- as.matrix(x)
    mod_x <- as.matrix(mod_x)
    y <- as.matrix(y)
    J <- dim(as.array(unique(y)))[1]
    if ( J <= 3 ){
        stop("This function is only available for models with more than 3 outcome
                variables.")
    }
    if (dim(y)[2] != 1){
        stop("parameter y should be a column vector")
    }
    if ( any(!all(y == floor(y)))){
        stop("each entry of y must be a integer")
    }
    if ( !all(is.numeric(x))){
        stop("each entry in x must be numeric")
    }
    if ( !all(is.numeric(mod_x))){
        stop("each entry in mod_x must be numeric")
    }
    if ( length(p) != 1){
        stop("parameter p must be scalar")
    }
    if (any(p < 0 | p > 1)){
        stop("parameter p must be between 0 to 1")
    }
    N <- dim(model$beta_draws)[2]
    m <- (N)/(1.25)
    burn <- 0.25 * m
    n <- dim(x)[1]
    k <- dim(x)[2]
    beta_burnt <- model$beta_draws[, (burn + 1):N]
    delta_burnt <- model$delta_draws[, (burn + 1):N]
    expdelta_burnt <- exp(delta_burnt)
    gammacp_cov <- array(0, dim = c(J-1, m))
    for (j in 2:(J-1)) {
            gammacp_cov[j, ] <- sum(expdelta_burnt[1:(j - 1), 1])
    }
    mu <- 0
    sigma <- 1
    old_prob <- array(0, dim = c(n, m, J))
    new_prob <- array(0, dim = c(n, m, J))
    old_comp <- array(0, dim = c(n, m, (J-1)))
    new_comp <- array(0, dim = c(n, m, (J-1)))
    for (j in 1:(J-1)) {
        for (i in 1:m) {
            for (k in 1:n) {
                old_comp[k, i, j] <- alcdf((gammacp_cov[j, i] - (x[k, ] %*% beta_burnt[, i])), mu, sigma, p)
                new_comp[k, i, j] <- alcdf((gammacp_cov[j, i] - (mod_x[k, ] %*% beta_burnt[, i])), mu, sigma, p)
            }
            if (j == 1) {
                old_prob[, i, j] <- old_comp[, i, j]
                new_prob[, i, j] <- new_comp[, i, j]
            }
            else {
                old_prob[, i, j] <- old_comp[, i, j] - old_prob[, i, (j-1)]
                new_prob[, i, j] <- new_comp[, i, j] - new_prob[, i, (j-1)]
            }
        }
    }
    old_prob[, , J] = 1 - old_comp[, , (J-1)]
    new_prob[, , J] = 1 - new_comp[, , (J-1)]
    diff_prob <- new_prob - old_prob
    avgdiff_prob = (colMeans(diff_prob, dims = 2))
    result <- list("avgdiff_prob" = avgdiff_prob)

    return(result)
}
