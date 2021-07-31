#' Bayesian Quantile Regression for Ordinal Models
#' with more than 3 outcomes
#'
#' This function estimates the Bayesian Quantile Regression for ordinal models with
#' more than 3 outcomes and reports the posterior mean and posterior standard deviations
#' of \eqn{(\beta, \delta)}.
#'
#' @param y         dependent variable i.e. ordinal outcome values.
#' @param x         covariate matrix of dimension \eqn{(n x k)} including a column of ones with or without column names.
#' @param mc        number of MCMC iterations, post burn-in.
#' @param p         quantile level or skewness parameter, p in (0,1).
#' @param tune      tuning parameter.
#'
#' @details
#' Function implements the Bayesian quantile regression for
#' ordinal models with more than 3 outcomes using a combination of Gibbs sampling
#' procedure and Metropolis-Hastings algorithm.
#'
#' Function initializes prior and then iteratively
#' samples \eqn{\beta}, \eqn{\delta} and latent variable z.
#' Burn-in is taken as \eqn{0.25*mc} and \eqn{N = burn}-\eqn{in + mc}.
#'
#' @return Returns a list with components:
#' \itemize{
#' \item{\code{post_mean_beta}: }{a vector with mean of sampled
#'  \eqn{\beta} for each covariate.}
#'  \item{\code{post_mean_beta}: }{a vector with mean of sampled
#'  \eqn{\beta} for each covariate.}
#'  \item{\code{post_mean_delta}: }{a vector with mean of sampled
#'  \eqn{\delta} for each cut point.}
#'  \item{\code{post_std_beta}: }{a vector with standard deviation
#'  of sampled \eqn{\beta} for each covariate.}
#'  \item{\code{post_std_delta}: }{a vector with standard deviation
#'  of sampled \eqn{\delta} for each cut point.}
#'  \item{\code{gamma}: }{a vector of cut points including Inf and
#'  -Inf.}
#'  \item{\code{catt}}
#'  \item{\code{acceptance_rate}: }{a scalar to judge the acceptance
#'  rate of samples.}
#'  \item{\code{DIC_result}: }{results of the DIC criteria.}
#'  \item{\code{beta_draws}: }{a matrix with all sampled values for \eqn{\beta}.}
#'  \item{\code{delta_draws}: }{a matrix with all sampled values for \eqn{\delta}.}
#' }
#'
#' @references Rahman, M. A. (2016). “Bayesian
#' Quantile Regression for Ordinal Models.”
#' Bayesian Analysis, 11(1): 1-24.
#'
#' Yu, K. and Moyeed, R. A. (2001). “Bayesian Quantile Regression.” Statistics and
#' Probability Letters, 54(4): 437–447.
#'
#' Casella, G., George E. I. (1992). “Explaining the Gibbs Sampler."
#' The American Statistician, 46(3): 167-174.
#'
#' Geman, S., and Geman, D. (1984). “Stochastic Relaxation,
#' Gibbs Distributions, and the Bayesian Restoration of Images."
#' IEEE Transactions an Pattern Analysis and Machine Intelligence,
#' 6(6): 721-741.
#'
#'  Chib, S., Greenberg E. (1995). "Understanding the Metropolis-Hastings
#'  Algorithm." The American Statistician, 49(4): 327-335.
#'
#'  Hastings, W.K. (1970). "Monte Carlo Sampling Methods Using
#'  Markov Chains and Their Applications." Biometrika, 57: 1317-1340.
#'
#' @importFrom "MASS" "ginv"
#' @importFrom "tcltk" "tkProgressBar" "setTkProgressBar"
#' @importFrom "stats" "sd"
#' @importFrom "stats" "rnorm" "qnorm"
#' @seealso tcltk, \link[stats]{rnorm}, \link[stats]{qnorm},
#' \link[MASS]{ginv}, Gibbs sampler
#' @examples
#'  set.seed(101)
#'  data("data25j4")
#'  x <- data25j4$x
#'  y <- data25j4$y
#'  p <- 0.25
#'  output <- quan_regg3(y, x, mc = 50, p, 0.89)
#'
#'
#'  # Number of burn-in draws: 12.5
#'  # Number of retained draws: 50
#'  # Summary of estimated beta:
#'
#'  #          Bayes Estimate Std Value Upper Credible Lower Credible
#'  # beta_0          0.7669    0.7669         1.3640         0.0862
#'  # beta_1          0.9110    0.9110         1.2914         0.4157
#'  # beta_2          3.2143    3.2143         4.0231         2.0865
#'  # delta_1         0.6665    0.2718         0.9781         0.2646
#'  # delta_2         0.8351    0.3719         1.1934         0.0040
#'
#'  # Model acceptance rate: 36
#'  # Model DIC: 732.6
#'
#' @export
quan_regg3 <- function(y, x, mc = 15000, p, tune = 0.1) {
    cols <- colnames(x)
    names(x) <- NULL
    names(y) <- NULL
    x <- as.matrix(x)
    y <- as.matrix(y)
    if (dim(y)[2] != 1){
        stop("parameter y should be a column vector")
    }
    if ( any(!all(y == floor(y)))){
        stop("each entry of y must be a integer")
    }
    if ( !all(is.numeric(x))){
        stop("each entry in x must be numeric")
    }
    if ( ( length(mc) != 1) || (length(tune) != 1 )){
        if (length(mc) != 1){
            stop("parameter mc must be scalar")
        }
        else{
            stop("parameter tune must be scalar")
        }
    }
    if (!is.numeric(mc)){
        stop("parameter mc must be a numeric")
    }
    if ( length(p) != 1){
        stop("parameter p must be scalar")
    }
    if (any(p < 0 | p > 1)){
        stop("parameter p must be between 0 to 1")
    }
    if (!is.numeric(tune)){
        stop("parameter tune must be numeric")
    }
    if (any(tune < 0 | tune > 1)){
        stop("parameter tune must be between 0 to 1")
    }
    J <- dim(as.array(unique(y)))[1]
    if ( J <= 3 ){
        warning("This function is for more than 3 outcome
                variables. We are switching to quan_reg3")
        ans <- quan_reg3(y, x, mc = mc, p)
    }
    n <- dim(x)[1]
    k <- dim(x)[2]
    burn <- 0.25 * mc
    yprob <- array(0, dim = c(n, J))
    for (i in 1:n) {
        yprob[i, y[i]] <- 1
    }
    yprob <- colSums(yprob) / n
    gam <- qnorm(cumsum(yprob[1:(J - 1)]))
    deltain <- t(log(gam[2:(J - 1)] - gam[1:(J - 2)]))

    b0 <- array(0, dim = c(k, 1))
    B0 <- diag(k)
    invb0 <- ginv(B0)
    invb0b0 <- invb0 %*% b0
    d0 <- array(0, dim = c(J-2, 1))
    D0 <- 0.25 * diag(J - 2)
    N <- burn + mc
    beta_draws <- array (0, dim = c(k, N))
    delta_draws <- array(0, dim = c(J - 2, N))
    ytemp <- y - 1.5
    beta_draws[, 1] <- mldivide( (t(x) %*% (x)), (t(x) %*% ytemp))
    delta_draws[, 1] <- deltain
    w <- array( (abs(rnorm(n, mean = 2, sd = 1))), dim = c (n, 1))
    z <- array( (rnorm(n, mean = 0, sd = 1)), dim = c(n, 1))

    theta <- (1 - (2 * p)) / (p * (1 - p))
    tau <- sqrt(2 / (p * (1 - p)))
    tau2 <- tau ^ 2
    lambda <- 0.5
    lengthp <- length(p)
    cri0     <- 1;
    cri1     <- 0.001;
    stepsize <- 1;
    maxiter  <- 10;
    h        <- 0.002;
    dh       <- 0.0002;
    sw       <- 20;
    minimize <- qrminfundtheorem(deltain, y, x,
                                 beta_draws[, 1], cri0, cri1,
                                 stepsize, maxiter, h, dh, sw, p)

    Dhat <- -ginv(minimize$H) * 3

    mhacc <- 0
    total <- mc
    pb <- tkProgressBar(title = "Simulation in Progress",
                        min = 0, max = N, width = 300)
    for (i in 2:N) {
        p1 <- drawbetag3(z, x, w, tau2, theta, invb0, invb0b0)
        beta_draws[, i] <- p1
        w <- drawwg3(z, x, beta_draws[, i], tau2, theta, lambda)
        delta0 <- delta_draws[, (i - 1)]

        deltarw <- drawdeltag3(y, x, beta_draws[, i], delta0, d0, D0, tune, Dhat, p)
        delta_draws[, i] <- deltarw$deltareturn
        z <- drawlatentg3(y, x, beta_draws[, (i - 1)], w, theta,
                          tau2, delta_draws[, i])

        if (i > burn) {
            mhacc <- mhacc + deltarw$accept
        }
        setTkProgressBar(pb, i,
                         label = paste(round( (i / N) * 100, 0), "% done"))
    }
    close(pb)
    post_mean_beta <- rowMeans(beta_draws[, (burn + 1):N])
    post_mean_delta <- rowMeans(delta_draws[, (burn + 1):N])
    post_std_beta <- apply(beta_draws[, (burn + 1):N], 1, sd)
    post_std_delta <- apply(delta_draws[, (burn + 1):N], 1, sd)
    gammacp <- array(0, dim = c(J - 1, 1))
    expdelta <- exp(post_mean_delta)
    for (j in 2:(J - 1)) {
        gammacp[j] <- sum(expdelta[1:(j - 1)])
    }

    acceptrate <- (mhacc / mc) * 100

    xbar <- colMeans(x)
    catt <- array(0, dim = c(J))
    catt[1] <- alcdfstdg3( (0 - (xbar %*% (post_mean_beta))), p)
    for (j in 2:(J - 1)) {
        catt[j] <- (alcdfstdg3(gammacp[j] - (xbar %*% (post_mean_beta)), p) -
                        alcdfstdg3(gammacp[(j - 1)] - (xbar %*% (post_mean_beta)),
                                   p))
    }
    catt[J] <- 1 - alcdfstdg3(gammacp[(J - 1)] - (xbar %*% (post_mean_beta)), p)
    DIC_result <- devianceg3(y, x, delta_draws, burn, N,
                             post_mean_beta, post_mean_delta, beta_draws, p)
    post_mean_beta <- array(post_mean_beta, dim = c(k, 1))
    post_std_beta <- array(post_mean_beta, dim = c(k, 1))
    post_mean_delta <- array(post_mean_delta, dim = c(J-2, 1))
    post_std_delta <- array(post_std_delta, dim = c(J-2, 1))
    upper_credible_beta <- array(apply(beta_draws[, (burn + 1):N], 1, quantile, c(0.975)), dim = c(k, 1))
    lower_credible_beta <- array(apply(beta_draws[, (burn + 1):N], 1, quantile, c(0.025)), dim = c(k, 1))
    upper_credible_delta <- array(apply(delta_draws[, (burn + 1):N], 1, quantile, c(0.975)), dim = c(J-2, 1))
    lower_credible_delta <- array(apply(delta_draws[, (burn + 1):N], 1, quantile, c(0.025)), dim = c(J-2, 1))
    beta_combine <- cbind(post_mean_beta, post_std_beta, upper_credible_beta, lower_credible_beta)
    delta_combine <- cbind(post_mean_delta, post_std_delta, upper_credible_delta, lower_credible_delta)
    summary <- rbind(beta_combine, delta_combine)
    name <- list('Bayes Estimate', 'Std Value', 'Upper Credible', 'Lower Credible')
    dimnames(summary)[[2]] <- name
    dimnames(summary)[[1]] <- letters[1:(k+J-2)]
    j <- 1
    if (is.null(cols)) {
        rownames(summary)[j] <- c('Intercept')
        for (i in paste0("beta_",1:k-1)) {
            rownames(summary)[j] = i
            j = j + 1
        }
    }
    else {
        for (i in cols) {
            rownames(summary)[j] = i
            j = j + 1
        }
    }
    for (i in paste0("delta_",1:(J-2))) {
        rownames(summary)[j] = i
        j = j + 1
    }
    print(noquote(paste0("Number of burn-in draws: ", burn)))
    print(noquote(paste0("Number of retained draws: ", mc)))
    print(noquote("Summary of estimated beta: "))
    cat("\n")
    print(round(summary, 4))
    cat("\n")
    print(noquote(paste0("Model acceptance rate: ", round(acceptrate, 2))))
    print(noquote(paste0("Model DIC: ", round(DIC_result$DIC, 2))))

    result <- list("post_mean_beta" = post_mean_beta,
                   "post_mean_delta" = post_mean_delta,
                   "post_std_beta" = post_std_beta,
                   "post_std_delta" = post_std_delta,
                   "gamma" = gammacp,
                   "catt" = catt,
                   "acceptancerate" = acceptrate,
                   "DIC_result" = DIC_result,
                   "beta_draws" = beta_draws,
                   "delta_draws" = delta_draws)
    return(result)
}
#' Minimize the negative of log-likelihood
#'
#' This function minimizes the negative of the log-likelihood for an
#' ordinal quantile model with respect to the cut-points \eqn{\delta} using the
#' Fundamental Theorem of Calculus.
#'
#' @param deltain   initialization of cut-points.
#' @param y         dependent variable i.e. ordinal outcome values.
#' @param x         covariate matrix of dimension \eqn{(n x k)} including a column of ones.
#' @param beta      column vector of coefficients of dimension \eqn{(k x 1)}.
#' @param cri0      initial criterion, \eqn{cri0 = 1}.
#' @param cri1      criterion lies between (0.001 to 0.0001).
#' @param stepsize  learning rate lies between (0.1, 1).
#' @param maxiter   maximum number of iteration.
#' @param h         change in value of each \eqn{\delta}, holding other \eqn{\delta}
#'                  constant for first derivatives.
#' @param dh        change in each value of \eqn{\delta}, holding other \eqn{\delta} constant
#'                  for second derivaties.
#' @param sw        iteration to switch from BHHH to inv(-H) algorithm.
#' @param p         quantile level or skewness parameter, p in (0,1).
#'
#' @details
#' First derivative from first principle
#' \deqn{dy/dx=[f(x+h)-f(x-h)]/2h}
#'
#' Second derivative from First principle
#'
#' \deqn{f'(x-h)=(f(x)-f(x-h))/h}
#'
#' \deqn{f''(x)= [{(f(x+h)-f(x))/h} - (f(x)-f(x-h))/h]/h}
#'
#'       \deqn{= [(f(x+h)+f(x-h)-2 f(x))]/h^2}
#'
#' cross partial derivatives
#'
#' \deqn{f(x) = [f(x+dh,y)-f(x-dh,y)]/2dh}
#'
#' \deqn{f(x,y)=[{(f(x+dh,y+dh) - f(x+dh,y-dh))/2dh} - {(f(x-dh,y+dh) -
#' f(x-dh,y-dh))/2dh}]/2dh}
#'
#' \deqn{= 0.25* [{(f(x+dh,y+dh)-f(x+dh,y-dh))} -{(f(x-dh,y+dh)-f(x-dh,y-dh))}]/dh2}
#'
#' @return Returns a list with components
#' \itemize{
#' \item{\code{dmin}: }{a vector with cutpoints that minimize the log-likelihood function.}
#' \item{\code{sumlogl}: }{a scalar with sum of log-likelihood values.}
#' \item{\code{logl}: }{a vector with log-likelihood values.}
#' \item{\code{G}: }{a gradient vector, \eqn{(n x k)} matrix with i-th row as the score
#' for the i-th unit.}
#' \item{\code{H}: }{represents Hessian matrix.}
#' }
#'
#' @importFrom "MASS" "ginv"
#' @importFrom "pracma" "mldivide"
#' @references Rahman, M. A. (2016). “Bayesian
#' Quantile Regression for Ordinal Models.”
#' Bayesian Analysis, 11(1): 1-24.
#'
#'
#' @seealso differential calculus, functional maximization,
#' \link[MASS]{ginv}, \link[pracma]{mldivide}
#' @examples
#' set.seed(101)
#' deltain <- c(-0.002570995,  1.044481071)
#' data("data25j4")
#' x <- data25j4$x
#' y <- data25j4$y
#' p <- 0.25
#' beta <- c(0.3990094, 0.8168991, 2.8034963)
#' cri0     <- 1
#' cri1     <- 0.001
#' stepsize <- 1
#' maxiter  <- 10
#' h        <- 0.002
#' dh       <- 0.0002
#' sw       <- 20
#' output <- qrminfundtheorem(deltain, y, x, beta, cri0, cri1, stepsize, maxiter, h, dh, sw, p)
#'
#' # deltamin
#' #   0.5225109 0.8782456
#' # negsum
#' #   276.1712
#' # logl
#' #   -1.84297177
#' #   -1.08096041
#' #   -1.52607866
#' #   -0.99466831
#' #   -1.00206354
#' #   -1.02615671
#' #   -1.86717309
#' #   -0.23775078
#' #   -0.96491711
#' #   -0.76172383.. soon
#' # G
#' #    0.803892784  0.00000000
#' #   -0.420190546  0.72908381
#' #   -0.421776117  0.72908341
#' #   -0.421776117 -0.60184063
#' #   -0.421776117 -0.60184063
#' #    0.151489598  0.86175120
#' #    0.296995920  0.96329114
#' #   -0.421776117  0.72908341
#' #   -0.340103190 -0.48530164
#' #    0.000000000  0.00000000
#' #   -0.421776117 -0.60184063.. soon
#' # H
#' #   -100.08669 -26.85551
#' #   -26.85551 -89.74705
#' # checkoutput
#' #   0    0    0    0    0    0    0 ... soon
#'
#' @export
qrminfundtheorem <- function(deltain, y, x, beta, cri0, cri1,
                             stepsize, maxiter, h, dh, sw, p) {
    if ( !all(is.numeric(deltain))){
        stop("each entry in deltain must be numeric")
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
    if ( !all(is.numeric(beta))){
        stop("each entry in beta must be numeric")
    }
    if (length(cri0) != 1){
        stop("parameter cri0 must be scalar")
    }
    if (length(cri1) != 1){
        stop("parameter cri1 must be scalar")
    }
    if (length(stepsize) != 1){
        stop("parameter stepsize must be scalar")
    }
    if (length(maxiter) != 1){
        stop("parameter maxiter must be scalar")
    }
    if (length(h) != 1){
        stop("parameter h must be scalar")
    }
    if (length(dh) != 1){
        stop("parameter dh must be scalar")
    }
    if (length(sw) != 1){
        stop("parameter sw must be scalar")
    }
    if (any(p < 0 | p > 1)){
        stop("parameter p must be between 0 to 1")
    }
    n <- length(y)
    d <- length(deltain)
    storevn <- array(0, dim = c(n, d))
    storevp <- array(0, dim = c(n, d))
    checkoutput <- array(0, dim = c(n, 3 * d + 1))

    cri <- cri0
    der <- array(0, dim = c(n, d))
    dh2 <- dh ^ 2

    jj <- 0
    while ( ( cri > cri1 ) && ( jj < maxiter )) {
        jj <- jj + 1
        number <- qrnegloglikensum(deltain, y, x, beta, p)
        vo <- -number$nlogl
        deltao <- deltain
        for (i in 1:d) {
            deltain[i] <- deltain[i] - h
            number1 <- qrnegloglikensum(deltain, y, x, beta, p)
            vn <- -number1$nlogl
            deltain <- deltao

            storevn[, i] <- vn
            deltain[i] <- deltain[i] + h
            number2 <- qrnegloglikensum(deltain, y, x, beta, p)
            vp <- -number2$nlogl
            deltain <- deltao

            storevp[, i] <- vp
            der[, i] <- ( (0.5 * (vp - vn))) / h;
        }
        hess <- array(0, dim = c(d, d))
        i <- 1
        j <- 1
        while (j <= d) {
            while (i <= j) {
                if (i == j) {
                    deltain[i] <- deltain[i] + dh
                    number3 <- qrnegloglikensum(deltain, y, x, beta, p)
                    vp2 <- -number3$nlogl
                    deltain <- deltao

                    deltain[i] <- deltain[i] - dh
                    number4 <- qrnegloglikensum(deltain, y, x, beta, p)
                    vn2 <- -number4$nlogl
                    deltain <- deltao

                    hess[i, j] <- sum( (vp2 + vn2 - (2 * vo)) / dh2)
                }
                else{
                    f <- c(i, j)
                    deltain[f] <- deltain[f] + dh
                    number5 <- qrnegloglikensum(deltain, y, x, beta, p)
                    vpp <- -number5$nlogl
                    deltain <- deltao

                    deltain[f] <- deltain[f] - dh
                    number6 <- qrnegloglikensum(deltain, y, x, beta, p)
                    vnn <- -number6$nlogl
                    deltain <- deltao

                    deltain[f] <- deltain[f] + c(dh, -dh)
                    number7 <- qrnegloglikensum(deltain, y, x, beta, p)
                    vpn <- -number7$nlogl
                    deltain <- deltao

                    deltain[f] <- deltain[f] + c(-dh, dh)
                    number8 <- qrnegloglikensum(deltain, y, x, beta, p)
                    vnp <- -number8$nlogl
                    deltain <- deltao

                    hess[i, j] <- 0.25 * sum( (vpp + vnn - vpn - vnp) / dh2)
                }
                i <- (i + 1)
            }
            i <- 1
            j <- (j + 1)
        }
        hess <- diag(1, nrow = d, ncol = d) * hess +
            (1 - diag(1, nrow = d, ncol = d)) * (hess + t(hess))
        cri <- sum(abs(colSums(der)))
        invder <- ginv(t(der) %*% (der))
        ddeltabhhh <- invder %*% ( (colSums(der)))
        ddeltahess <- mldivide(-hess, (colSums(der)))

        ddelta <- ( ( (1 - min(1, max(0, jj - sw))) * ddeltabhhh) +
                        (min(1, max(0, jj - sw)) * ddeltahess))
        deltain <- deltain + stepsize * t(ddelta)

        if (jj == maxiter){
        }
    }
    deltamin <- deltain
    number9 <- qrnegloglikensum(deltamin, y, x, beta, p)
    logl <- -number9$nlogl
    negsum <- number9$negsumlogl
    G <- der
    H <- hess
    rt <- list("deltamin" = deltamin,
               "negsum" = negsum,
               "logl" = logl,
               "G" = G,
               "H" = H,
               "checkoutput" = checkoutput)
    return(rt)
}

#' Negative log-likelihood for Ordinal Models with more than 3 outcomes
#'
#' Function for calculating negative log-likelihood for Ordinal models with
#' more than 3 outcomes.
#'
#' @param deltain   initialization of cut-points.
#' @param y         dependent variable i.e. ordinal outcome values.
#' @param x         covariate matrix of dimension \eqn{(n x k)} including a column of ones.
#' @param beta      column vector of coefficients of dimension \eqn{(k x 1)}.
#' @param p         quantile level or skewness parameter, p in (0,1).
#'
#' @details
#' Computes the negative of the log-likelihood function using the
#' asymmetric Laplace distribution over the iid random
#' variables.
#'
#' @return Returns a list with components
#' \itemize{
#' \item{\code{nlogl}: }{a vector with likelihood values.}
#' \item{\code{negsumlogl}: }{a scalar with value of negative log-likelihood.}
#' }
#'
#' @references Rahman, M. A. (2016). “Bayesian
#' Quantile Regression for Ordinal Models.”
#' Bayesian Analysis, 11(1): 1-24.
#'
#'
#' @seealso likelihood maximization
#' @examples
#' set.seed(101)
#' deltain <- c(-0.002570995, 1.044481071)
#' data("data25j4")
#' x <- data25j4$x
#' y <- data25j4$y
#' p <- 0.25
#' beta <- c(0.3990094, 0.8168991, 2.8034963)
#' output <- qrnegloglikensum(deltain, y, x, beta, p)
#'
#' # nlogl
#' #   2.28675506
#' #   0.88417292
#' #   1.23613687
#' #   0.93127471
#' #   0.93866993
#' #   0.97694495
#' #   1.04934754
#' #   1.57723130
#' #   0.19231567
#' #   0.96491711... soon
#'
#' # negsumlogl
#' #   285.7158
#'
#' @export
qrnegloglikensum <- function(deltain, y, x, beta, p) {
    if ( !all(is.numeric(deltain))){
        stop("each entry in deltain must be numeric")
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
    if ( !all(is.numeric(beta))){
        stop("each entry in beta must be numeric")
    }
    if (any(p < 0 | p > 1)){
        stop("parameter p must be between 0 to 1")
    }
    J <- dim(as.array(unique(y)))[1]
    n <- dim(x)[1]
    lnpdf <- array(0, dim = c(n, 1))
    expdelta <- exp(deltain)
    q <- (dim(expdelta)[1]) + 1
    gammacp <- array(0, dim = c(q, 1))

    for (j in 2:(J - 1)) {
        gammacp[j] <- sum(expdelta[1:(j - 1)])
    }
    allgammacp <- t(c(-Inf, gammacp, Inf))
    mu <- x %*% beta
    for (i in 1:n) {
        meanp <- mu[i]
        if (y[i] == 1) {
            lnpdf[i] <- log(alcdf(0, meanp, 1, p))
        }
        else if (y[i] == J) {
            lnpdf[i] <- log(1 - alcdf(allgammacp[J], meanp, 1, p))
        }
        else{
            lnpdf[i] <- log( (alcdf(allgammacp[y[i] + 1], meanp, 1, p) -
                                  alcdf(allgammacp[y[i]], meanp, 1, p)))
        }
    }
    nlogl <- -lnpdf
    negsumlogl <- sum(nlogl)
    respon <- list("nlogl" = nlogl,
                   "negsumlogl" = negsumlogl)
    return(respon)
}
#' Samples \eqn{\beta} for an Ordinal Model
#' with more than 3 outcomes
#'
#' This function samples \eqn{\beta} from its conditional
#' posterior distribution for an ordinal model with more than 3
#' outcomes.
#'
#' @param z         Gibbs draw of latent response variable, a column vector.
#' @param x         covariate matrix of dimension \eqn{(n x k)} including a column of ones.
#' @param w         latent weights, row vector.
#' @param tau2      2/(p(1-p)).
#' @param theta     (1-2p)/(p(1-p)).
#' @param invB0     inverse of prior covariance matrix of normal distribution.
#' @param invB0b0   prior mean pre-multiplied by invB0.
#'
#' @details
#' Function samples a vector of \eqn{\beta} from a multivariate
#' normal distribution.
#'
#' @return Returns a column vector of \eqn{\beta}
#' from a multivariate normal distribution.
#'
#' @references Rahman, M. A. (2016). “Bayesian
#' Quantile Regression for Ordinal Models.”
#' Bayesian Analysis, 11(1): 1-24.
#'
#' Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988).
#' “The New S Language. Wadsworth & Brooks/Cole.”
#'
#' Casella, G., George E. I. (1992). “Explaining the Gibbs Sampler."
#' The American Statistician, 46(3): 167-174.
#'
#' Geman, S., and Geman, D. (1984). “Stochastic Relaxation,
#' Gibbs Distributions, and the Bayesian Restoration of Images."
#' IEEE Transactions an Pattern Analysis and Machine Intelligence,
#' 6(6): 721-741.
#'
#' @importFrom "MASS" "mvrnorm" "ginv"
#' @seealso Gibbs sampling, normal distribution,
#' \link[MASS]{ginv},  \link[MASS]{mvrnorm}
#' @examples
#' set.seed(101)
#' data("data25j4")
#' x <- data25j4$x
#' p <- 0.25
#' n <- dim(x)[1]
#' k <- dim(x)[2]
#' w <- array( (abs(rnorm(n, mean = 2, sd = 1))), dim = c (n, 1))
#' theta <- 2.666667
#' tau2 <- 10.66667
#' z <- array( (rnorm(n, mean = 0, sd = 1)), dim = c(n, 1))
#' b0 <- array(0, dim = c(k, 1))
#' B0 <- diag(k)
#' invB0 <- matrix(c(
#'      1, 0, 0,
#'      0, 1, 0,
#'      0, 0, 1),
#'      nrow = 3, ncol = 3, byrow = TRUE)
#' invB0b0 <- invB0 %*% b0
#' output <- drawbetag3(z, x, w, tau2, theta, invB0, invB0b0)
#'
#' # output
#' #   -5.00085 -2.500239 -0.5525601
#' @export
drawbetag3 <- function(z, x, w, tau2, theta, invB0, invB0b0) {
    if ( !all(is.numeric(z))){
        stop("each entry in z must be numeric")
    }
    if ( !all(is.numeric(x))){
        stop("each entry in x must be numeric")
    }
    if ( !all(is.numeric(w))){
        stop("each entry in w must be numeric")
    }
    if ( length(tau2) != 1){
        stop("parameter tau2 must be scalar")
    }
    if ( !all(is.numeric(tau2))){
        stop("each entry in tau2 must be numeric")
    }
    if ( length(theta) != 1){
        stop("parameter theta must be scalar")
    }
    if ( !all(is.numeric(theta))){
        stop("each entry in theta must be numeric")
    }
    if ( !all(is.numeric(invB0))){
        stop("each entry in invB0 must be numeric")
    }
    if ( !all(is.numeric(invB0b0))){
        stop("each entry in invB0b0 must be numeric")
    }
    n <- dim(x)[1]
    k <- dim(x)[2]
    meancomp <- array(0, dim = c(n, k))
    varcomp <- array(0, dim = c(k, k, n))
    q <- array(0, dim = c(1, k))
    eye <- diag(k)
    for (j in 1:n){
        meancomp[j, ] <- (x[j, ] * (z[j] - (theta * w[j]))) / (tau2 * w[j])
        varcomp[, , j] <- ( (x[j, ]) %*% (t(x[j, ]))) / (tau2 * w[j])
    }
    Btilde <- ginv(invB0 + rowSums(varcomp, dims = 2))
    btilde <- Btilde %*% (invB0b0 + colSums(meancomp))
    L <- t(chol(Btilde))
    beta <- t(btilde) + t(L %*% ( (mvrnorm(n = 1, mu = q, Sigma = eye))))
    return(beta)
}
#' Samples the latent weight w for an Ordinal Model
#' with more than 3 outcomes
#'
#' This function samples the latent weight w from a Generalized
#' inverse-Gaussian distribution (GIG) for an ordinal model with more
#' than 3 outcomes.
#'
#' @param z         Gibbs draw of latent response variable, a column vector.
#' @param x         covariate matrix of dimension \eqn{(n x k)} including a column of ones.
#' @param beta      Gibbs draw of coefficients of dimension \eqn{(k x 1)}.
#' @param tau2      2/(p(1-p)).
#' @param theta     (1-2p)/(p(1-p)).
#' @param lambda    index parameter of GIG distribution which is equal to 0.5
#'
#' @details
#' Function samples a vector of the latent weight w from a GIG distribution.
#'
#' @return Returns a column vector of w from a
#' GIG distribution.
#'
#' @references Rahman, M. A. (2016). “Bayesian
#' Quantile Regression for Ordinal Models.”
#' Bayesian Analysis, 11(1): 1-24.
#'
#' Casella, G., George E. I. (1992). “Explaining the Gibbs Sampler."
#' The American Statistician, 46(3): 167-174.
#'
#' Geman, S., and Geman, D. (1984). “Stochastic Relaxation,
#' Gibbs Distributions, and the Bayesian Restoration of Images."
#' IEEE Transactions an Pattern Analysis and Machine Intelligence,
#' 6(6): 721-741.
#'
#' @importFrom "GIGrvg" "rgig"
#' @seealso GIGrvg, Gibbs sampling, \link[GIGrvg]{rgig}
#' @examples
#' set.seed(101)
#' z <- c(0.9812363, -1.09788, -0.9650175, 8.396556,
#'  1.39465, -0.8711435, -0.5836833, -2.792464,
#'  0.1540086, -2.590724, 0.06169976, -1.823058,
#'  0.06559151, 0.1612763, 0.161311, 4.908488,
#'  0.6512113, 0.1560708, -0.883636, -0.5531435)
#' x <- matrix(c(
#'      1, 1.4747905363, 0.167095186,
#'      1, -0.3817326861, 0.041879526,
#'      1, -0.1723095575, -1.414863777,
#'      1, 0.8266428137, 0.399722073,
#'      1, 0.0514888733, -0.105132425,
#'      1, -0.3159992662, -0.902003846,
#'      1, -0.4490888878, -0.070475600,
#'      1, -0.3671705251, -0.633396477,
#'      1, 1.7655601639, -0.702621934,
#'      1, -2.4543678120, -0.524068780,
#'      1,  0.3625025618,  0.698377504,
#'      1, -1.0339179063,  0.155746376,
#'      1,  1.2927374692, -0.155186911,
#'      1, -0.9125108094, -0.030513775,
#'      1,  0.8761233001,  0.988171587,
#'      1,  1.7379728231,  1.180760114,
#'      1,  0.7820635770, -0.338141095,
#'      1, -1.0212853209, -0.113765067,
#'      1,  0.6311364051, -0.061883874,
#'      1,  0.6756039688,  0.664490143),
#'      nrow = 20, ncol = 3, byrow = TRUE)
#' beta <- c(-1.583533, 1.407158, 2.259338)
#' tau2 <- 10.66667
#' theta <- 2.666667
#' lambda <- 0.5
#' output <- drawwg3(z, x, beta, tau2, theta, lambda)
#'
#' # output
#' #   0.16135732
#' #   0.39333080
#' #   0.80187227
#' #   2.27442898
#' #   0.90358310
#' #   0.99886987
#' #   0.41515947 ... soon
#'
#' @export
drawwg3 <- function(z, x, beta, tau2, theta, lambda) {
    if ( !all(is.numeric(z))){
        stop("each entry in z must be numeric")
    }
    if ( !all(is.numeric(x))){
        stop("each entry in x must be numeric")
    }
    if ( !all(is.numeric(beta))){
        stop("each entry in beta must be numeric")
    }
    if ( length(tau2) != 1){
        stop("parameter tau2 must be scalar")
    }
    if ( !all(is.numeric(tau2))){
        stop("each entry in tau2 must be numeric")
    }
    if ( length(theta) != 1){
        stop("parameter theta must be scalar")
    }
    if ( !all(is.numeric(theta))){
        stop("each entry in theta must be numeric")
    }
    if ( length(lambda) != 1){
        stop("parameter lambda must be scalar")
    }
    if ( !all(is.numeric(lambda))){
        stop("each entry in lambda must be numeric")
    }
    n <- dim(x)[1]
    tildeeta2 <- ( (theta ^ 2) / (tau2)) + 2
    tildelambda2 <- array(0, dim = c(n, 1))
    w <- array(0, dim = c(n, 1))
    for (j in 1:n) {
        tildelambda2[j, 1] <- ( (z[j] - (x[j, ] %*% beta)) ^ 2) / (tau2)
        if(is.na(tildelambda2[j, 1]) == TRUE){
            tildelambda2[j, 1] <- 0.5
        }
        w[j, 1] <- rgig(lambda = lambda,
                        chi = tildelambda2[j, 1],
                        psi = tildeeta2,
                        n = 1)
    }
    return(w)
}
#' Samples the Latent Variable z for an Ordinal Models
#' with more than 3 outcomes
#'
#' This function samples the latent variable z from a truncated
#' normal distribution for an ordinal model with more than 3 outcomes.
#'
#' @param y         dependent variable i.e. ordinal outcome values.
#' @param x         covariate matrix of dimension \eqn{(n x k)} including a column of ones.
#' @param beta      Gibbs draw of coefficients of dimension \eqn{(k x 1)}.
#' @param w         latent weights vector.
#' @param theta     (1-2p)/(p(1-p)).
#' @param tau2      2/(p(1-p)).
#' @param delta     row vector of cutpoints including -Inf and Inf.
#'
#' @details
#' Function samples the latent variable z from a truncated normal
#' distribution.
#'
#' @return Returns a column vector of values for latent variable, z.
#'
#' @references Albert, J. and Chib, S. (1993). “Bayesian Analysis of Binary and Polychotomous
#' Response Data.” Journal of the American Statistical Association, 88(422): 669–679.
#'
#' Casella, G., George E. I. (1992). “Explaining the Gibbs Sampler."
#' The American Statistician, 46(3): 167-174.
#'
#' Geman, S., and Geman, D. (1984). “Stochastic Relaxation,
#' Gibbs Distributions, and the Bayesian Restoration of Images."
#' IEEE Transactions an Pattern Analysis and Machine Intelligence,
#' 6(6): 721-741.
#'
#' @importFrom "truncnorm" "rtruncnorm"
#' @seealso Gibbs sampling, truncated normal distribution,
#' \link[truncnorm]{rtruncnorm}
#' @examples
#' set.seed(101)
#' data("data25j4")
#' x <- data25j4$x
#' y <- data25j4$y
#' p <- 0.25
#' beta <- c(0.3990094, 0.8168991, 2.8034963)
#' w <- 1.114347
#' theta <- 2.666667
#' tau2 <- 10.66667
#' delta <- c(-0.002570995,  1.044481071)
#' output <- drawlatentg3(y, x, beta, w, theta, tau2, delta)
#'
#' # output
#' #   0.6261896 3.129285 2.659578 8.680291
#' #   13.22584 2.545938 1.507739 2.167358
#' #   15.03059 -3.963201 9.237466 -1.813652
#' #   2.718623 -3.515609 8.352259 -0.3880043
#' #   -0.8917078 12.81702 -0.2009296 1.069133 ... soon
#'
#' @export
drawlatentg3 <- function(y, x, beta, w, theta, tau2, delta) {
    if (dim(y)[2] != 1){
        stop("parameter y should be a column vector")
    }
    if ( any(!all(y == floor(y)))){
        stop("each entry of y must be a integer")
    }
    if ( !all(is.numeric(x))){
        stop("each entry in x must be numeric")
    }
    if ( !all(is.numeric(beta))){
        stop("each entry in beta must be numeric")
    }
    if ( !all(is.numeric(w))){
        stop("each entry in w must be numeric")
    }
    if ( length(tau2) != 1){
        stop("parameter tau2 must be scalar")
    }
    if ( !all(is.numeric(tau2))){
        stop("each entry in tau2 must be numeric")
    }
    if ( length(theta) != 1){
        stop("parameter theta must be scalar")
    }
    if ( !all(is.numeric(theta))){
        stop("each entry in theta must be numeric")
    }
    if ( !all(is.numeric(delta))){
        stop("each entry in delta must be numeric")
    }
    J <- dim(as.array(unique(y)))[1]
    n <- dim(x)[1]
    z <- array(0, dim = c(1, n))
    expdelta <- exp(delta)
    q <- dim(expdelta)[1] + 1
    gammacp <- array(0, dim = c(q, 1))

    for (j in 2:(J - 1)) {
        gammacp[j] <- sum(expdelta[1:(j - 1)])
    }
    gammacp <- t(c(-Inf, gammacp, Inf))
    for (i in 1:n) {
        meanp <- (x[i, ] %*% (beta)) + (theta * w[i])
        std <- sqrt(tau2 * w[i])
        temp <- y[i]
        a <- gammacp[temp]
        b <- gammacp[temp + 1]
        z[1, i] <- rtruncnorm(n = 1, a = a, b = b,
                                  mean = meanp, sd = std)
        if(is.na(z[1, i]) == TRUE) {
            z[1, i] <- b
        }
    }
    return(z)
}
#' Samples the \eqn{\delta} for an Ordinal Model
#' with more than 3 outcomes
#'
#' This function samples the \eqn{\delta} using a
#' random-walk Metropolis-Hastings algorithm for an ordinal
#' model with more than 3 outcomes.
#'
#' @param y         dependent variable i.e. ordinal outcome values..
#' @param x         covariate matrix of dimension \eqn{(n x k)} including a column of ones.
#' @param beta      Gibbs draw of coefficients of dimension \eqn{(k x 1)}.
#' @param delta0    initial value for \eqn{\delta}.
#' @param d0        prior mean of normal distribution.
#' @param D0        prior variance-covariance matrix of normal distribution.
#' @param tune      tuning parameter.
#' @param Dhat      negative inverse Hessian from maximization of log-likelihood.
#' @param p         quantile level or skewness parameter, p in (0,1).
#'
#' @details
#' Samples the \eqn{\delta} using a random-walk Metropolis-Hastings algorithm.
#'
#' @return Returns a list with components
#' \itemize{
#'  \item{\code{deltaReturn}: }{a vector with \eqn{\delta} values using MH algorithm.}
#'   \item{\code{accept}: }{an indicator for acceptance of proposed value of \eqn{\delta}.}
#' }
#'
#' @references Rahman, M. A. (2016). “Bayesian
#' Quantile Regression for Ordinal Models.”
#' Bayesian Analysis, 11(1): 1-24.
#'
#'  Chib, S., Greenberg E. (1995). "Understanding the Metropolis-Hastings
#'  Algorithm." The American Statistician, 49(4): 327-335.
#'
#'  Hastings, W.K. (1970). "Monte Carlo Sampling Methods Using
#'  Markov Chains and Their Applications." Biometrika, 57: 1317-1340.
#'
#' @importFrom "stats" "rnorm"
#' @importFrom "pracma" "rand"
#' @importFrom "NPflow" "mvnpdf"
#' @seealso NPflow, Gibbs sampling, \link[NPflow]{mvnpdf}
#' @examples
#' set.seed(101)
#' data("data25j4")
#' x <- data25j4$x
#' y <- data25j4$y
#' p <- 0.25
#' beta <- c(0.3990094, 0.8168991, 2.8034963)
#' delta0 <- c(-0.9026915, -2.2488833)
#' d0 <- matrix(c(0, 0),
#'                  nrow = 2, ncol = 1, byrow = TRUE)
#' D0 <- matrix(c(0.25, 0.00, 0.00, 0.25),
#'                     nrow = 2, ncol = 2, byrow = TRUE)
#' tune <- 0.1
#' Dhat <- matrix(c(0.046612180, -0.001954257, -0.001954257, 0.083066204),
#'              nrow = 2, ncol = 2, byrow = TRUE)
#' p <- 0.25
#' output <- drawdeltag3(y, x, beta, delta0, d0, D0, tune, Dhat, p)
#'
#' # deltareturn
#' #   -0.9025802 -2.229514
#' # accept
#' #   1
#'
#' @export
drawdeltag3 <- function(y, x, beta, delta0, d0, D0, tune, Dhat, p){
    if (dim(y)[2] != 1){
        stop("parameter y should be a column vector")
    }
    if ( any(!all(y == floor(y)))){
        stop("each entry of y must be a integer")
    }
    if ( !all(is.numeric(x))){
        stop("each entry in x must be numeric")
    }
    if ( !all(is.numeric(beta))){
        stop("each entry in beta must be numeric")
    }
    if ( !all(is.numeric(delta0))){
        stop("each entry in delta0 must be numeric")
    }
    if ( !all(is.numeric(d0))){
        stop("each entry in d0 must be numeric")
    }
    if ( !all(is.numeric(D0))){
        stop("each entry in D0 must be numeric")
    }
    if ( length(tune) != 1){
        stop("parameter tune must be a scalar")
    }
    if (!is.numeric(tune)){
        stop("parameter tune must be numeric")
    }
    if (any(tune < 0 | tune > 1)){
        stop("parameter tune must be between 0 to 1")
    }
    if ( !all(is.numeric(Dhat))){
        stop("each entry in Dhat must be numeric")
    }
    if (any(p < 0 | p > 1)){
        stop("parameter p must be between 0 to 1")
    }
    J <- dim(as.array(unique(y)))[1]
    k <- (J - 2)
    L <- t(chol(Dhat))
    delta1 <- delta0 + tune * t(L %*% ( (rnorm(n = k, mean = 0, sd = 1))))
    num <-  qrnegloglikensum(delta1, y, x, beta, p)
    den <-  qrnegloglikensum(delta0, y, x, beta, p)
    pnum <- - (num$negsumlogl +
                   mvnpdf(x = matrix(delta1),
                          mean  = matrix(d0),
                          varcovM = D0,
                          Log = TRUE))
    pden <- - (den$negsumlogl +
                   mvnpdf(x = matrix(delta0),
                          mean  = matrix(d0),
                          varcovM = D0,
                          Log = TRUE))
    if((is.infinite(pnum) == TRUE) & (is.infinite(pden) == TRUE) & (pden < 0) & (pnum < 0)){
        pnum = -1000000
        pden = -1000000
    }
    if (log(rand(n = 1)) <= (pnum - pden)) {
        deltareturn <- delta1
        accept <- 1
    }
    else{
        deltareturn <- delta0
        accept <- 0
    }
    resp <- list("deltareturn" = deltareturn,
                 "accept" = accept)
    return(resp)
}
#' Deviance Information Criteria for Ordinal Models
#' with more than 3 outcomes
#'
#' Function for computing the Deviance Information Criteria for ordinal
#' models with more than 3 outcomes.
#'
#' @param y                dependent variable i.e. ordinal outcome values.
#' @param x                covariate matrix of dimension \eqn{(n x k)} including a column of ones.
#' @param post_mean_beta   mean value of \eqn{\beta} obtained from MCMC draws.
#' @param post_mean_delta  mean value of \eqn{\delta} obtained from MCMC draws.
#' @param beta_draws       MCMC draw of coefficients, dimension is \eqn{(k x N)}.
#' @param deltastore       MCMC draws of \eqn{\delta}.
#' @param p                quantile level or skewness parameter, p in (0,1).
#' @param burn             number of discarded MCMC iterations.
#' @param N                total number of samples, including the burn-in.
#'
#' @details
#' The Deviance is -2*(log likelihood) and has an important role in
#' statistical model comparision because of its relation with Kullback-Leibler
#' information criteria.
#'
#' @return Returns a list with components
#' \deqn{DIC = 2*avgdeviance - devpostmean}
#' \deqn{pd = avgdeviance - devpostmean}
#' \deqn{devpostmean = -2*(logLikelihood)}.
#'
#' @references Rahman, M. A. (2016). “Bayesian
#' Quantile Regression for Ordinal Models.”
#' Bayesian Analysis, 11(1): 1-24.
#'
#' Spiegelhalter, D. J., Best, N. G., Carlin B. P. and Linde A. (2002).
#' “Bayesian Measures of Model Complexity and Fit.” Journal of the
#' Royal Statistical Society B, Part 4: 583-639.
#'
#' Gelman, A., Carlin, J. B., Stern, H. S., and Rubin, D. B.
#' “Bayesian Data Analysis.” 2nd Edition, Chapman and Hall.
#'
#' @seealso  decision criteria
#' @examples
#' set.seed(101)
#' data("data25j4")
#' x <- data25j4$x
#' y <- data25j4$y
#' p <- 0.25
#' output <- quan_regg3(y, x, mc = 50, p, 0.75)
#' mc <- 50
#' deltastore <- output$delta_draws
#' burn <- 0.25*mc
#' N <- burn + mc
#' post_mean_beta <- output$post_mean_beta
#' post_mean_delta <- output$post_mean_delta
#' beta_draws <- output$beta_draws
#' deviance <- devianceg3(y, x, deltastore, burn, N,
#' post_mean_beta, post_mean_delta, beta_draws, p)
#'
#' # DIC
#' #   727.9137
#' # pd
#' #   78.241
#' # devpostmean
#' #   571.4317
#'
#' @export
devianceg3 <- function(y, x, deltastore, burn,
                       N, post_mean_beta, post_mean_delta, beta_draws, p){
    cols <- colnames(x)
    names(x) <- NULL
    names(y) <- NULL
    x <- as.matrix(x)
    y <- as.matrix(y)
    if (dim(y)[2] != 1){
        stop("parameter y should be a column vector")
    }
    if ( any(!all(y == floor(y)))){
        stop("each entry of y must be a integer")
    }
    if ( !all(is.numeric(x))){
        stop("each entry in x must be numeric")
    }
    if ( !all(is.numeric(deltastore))){
        stop("each entry in deltastore must be numeric")
    }
    if ( length(burn) != 1){
        stop("parameter burn must be scalar")
    }
    if ( length(N) != 1){
        stop("parameter N must be scalar")
    }
    if ( !all(is.numeric(post_mean_beta))){
        stop("each entry in post_mean_beta must be numeric")
    }
    if ( !all(is.numeric(post_mean_delta))){
        stop("each entry in post_mean_delta must be numeric")
    }
    if ( !all(is.numeric(beta_draws))){
        stop("each entry in beta must be numeric")
    }
    if (any(p < 0 | p > 1)){
        stop("parameter p must be between 0 to 1")
    }
    delta <- deltastore
    devpostmean <- array(0, dim = c(1))
    DIC <- array(0, dim = c(1))
    pd <- array(0, dim = c(1))
    ans <- qrnegloglikensum(post_mean_delta, y, x,
                            post_mean_beta, p)
    devpostmean <- 2 * ans$negsumlogl
    nsim <- dim(beta_draws[, (burn + 1):N])[1]
    devianceg3 <- array(0, dim = c(nsim, 1))
    for (j in 1:nsim) {
        temp <- qrnegloglikensum(delta[, (burn + j)],
                                 y, x, beta_draws[, (burn + j)], p)
        devianceg3[j, 1] <- 2 * temp$negsumlogl
    }
    avgdeviance <- mean(devianceg3)
    DIC <- 2 * avgdeviance - devpostmean
    pd <- avgdeviance - devpostmean
    result <- list("DIC" = DIC,
                   "pd" = pd,
                   "devpostmean" = devpostmean)
    return(result)
}
#' CDF of a standard Asymmetric Laplace Distribution
#'
#' This function computes the CDF of a standard asymmetric
#' Laplace distribution i.e. AL\eqn{(0, 1 ,p)}.
#'
#' @param x     scalar value.
#' @param p     quantile level or skewness parameter, p in (0,1).
#'
#' @details
#' Computes the CDF of a standard asymmetric Laplace distribution.
#' \deqn{CDF(x) = F(x) = P(X \le x)} where X is a
#' random variable that follows AL\eqn{(0, 1 ,p)}.
#'
#' @return Returns the probability value from the CDF of an asymmetric
#' Laplace distribution.
#'
#' @references Rahman, M. A. (2016). “Bayesian
#' Quantile Regression for Ordinal Models.”
#' Bayesian Analysis, 11(1): 1-24.
#'
#'  Koenker, R. and Machado, J. (1999). “Goodness  of  Fit  and  Related
#'  Inference  Processes  for  Quantile Regression.”
#'  Journal of American Statistics Association, 94(3): 1296-1309.
#'
#'  Keming, Y. and Zhang, J. (2005). “A Three-Parameter Asymmetric
#'  Laplace Distribution.” Communications in Statistics - Theory and Methods, 34(9):
#'  1867-1879.
#'
#' @seealso asymmetric Laplace distribution
#' @examples
#' set.seed(101)
#' x <-  -0.5428573
#' p <- 0.25
#' output <- alcdfstdg3(x, p)
#'
#' # output
#' #   0.1663873
#'
#' @export
alcdfstdg3 <- function(x, p) {
    if ( any(p < 0 | p > 1)){
        stop("parameter p must be between 0 to 1")
    }
    if (length(x) != 1){
        stop("parameter x must be scalar")
    }
    if (x <= 0) {
        z <- p * (exp( (1 - p) * x))
    }
    else {
        z <- 1 - (1 - p) * exp(-p * x )
    }
    return(z)
}
#' Asymmetric Laplace Distribution
#'
#' This function computes the cumulative distribution (CDF) for
#' an asymmetric Laplace distribution.
#'
#' @param x     scalar value.
#' @param mu    location parameter of ALD.
#' @param sigma scale parameter of ALD.
#' @param p     quantile or skewness parameter, p in (0,1).
#'
#' @details
#' Computes the cumulative distribution function for
#' the asymmetric Laplace distribution.
#' \deqn{CDF(x) = F(x) = P(X \le x)} where X is a
#' random variable
#'
#' @return Returns a scalar with cumulative probability value at
#' point 'x'.
#'
#' @references Rahman, M. A. (2016). “Bayesian
#' Quantile Regression for Ordinal Models.”
#' Bayesian Analysis, 11(1): 1-24.
#'
#'  Koenker, R. and Machado, J. (1999). “Goodness  of  Fit  and  Related
#'  Inference  Processes  for  Quantile Regression.”
#'  Journal of American Statistics Association, 94(3): 1296-1309.
#'
#'  Keming, Y. and Zhang, J. (2005). “A Three-Parameter Asymmetric
#'  Laplace Distribution.” Communications in Statistics - Theory and Methods, 34(9):
#'  1867-1879.
#'
#' @seealso cumulative distribution function, asymmetric Laplace distribution
#' @examples
#' set.seed(101)
#' x <- -0.5428573
#' mu <- 0.5
#' sigma <- 1
#' p <- 0.25
#' output <- alcdf(x, mu, sigma, p)
#'
#' # output
#' #   0.1143562
#'
#' @export
alcdf <- function(x, mu, sigma, p){
    if ( any(p < 0 | p > 1)){
        stop("parameter p must be between 0 to 1")
    }
    if ( ( length(mu) != 1) || (length(x) != 1 )
         || (length(sigma) != 1)){
        if (length(mu) != 1){
            stop("parameter mu must be scalar")
        }
        else if (length(x) != 1){
            stop("parameter x must be scalar")
        }
        else{
            stop("parameter sigma must be scalar")
        }
    }
    if (x <= mu) {
        z <- p * exp( (1 - p) * ( (x - mu) / sigma))
    }
    else {
        z <- 1 - ( (1 - p) * exp(-p * ( (x - mu) / sigma)))
    }
    return(z)
}
#' Trace Plots for Ordinal Models
#' with more than 3 outcomes
#'
#' This function generates trace plots of
#' MCMC samples for \eqn{(\beta ,\delta)} in the quantile
#' regression model with more than 3 outcomes.
#'
#' @param beta_draws      Gibbs draw of \eqn{\beta} vector of dimension \eqn{(k x N)}.
#' @param delta_draws     Gibbs draw of \eqn{\delta}.
#'
#' @details
#' Trace plot is a visual depiction of the values generated from the Markov chain
#' versus the iteration number.
#'
#' @return Returns trace plots for each element of \eqn{\beta}
#' and \eqn{\delta}.
#'
#' @references Rahman, M. A. (2016). “Bayesian
#' Quantile Regression for Ordinal Models.”
#' Bayesian Analysis, 11(1): 1-24.
#'
#' @importFrom "graphics" "plot"
#' @seealso traces in MCMC simulations
#' @examples
#' set.seed(101)
#' data("data25j4")
#' x <- data25j4$x
#' y <- data25j4$y
#' p <- 0.25
#' output <- quan_regg3(y, x, mc = 50, p, 0.75)
#' beta_draws <- output$beta_draws
#' delta_draws <- output$delta_draws
#' trace_plotg3(beta_draws, delta_draws)
#'
#' @export
trace_plotg3 <- function(beta_draws, delta_draws) {
    if ( !all(is.numeric(beta_draws))){
        stop("each entry in beta_draws must be numeric")
    }
    N <- dim(beta_draws)[2]
    k <- dim(beta_draws)[1]
    leng <- seq(0, N - 1)
    for (t in 1:k) {
        plot(beta_draws[t, ], xlab = '', ylab = '', type = "l", cex.main = 1.5, main = expression(paste(Trace~plot~beta)))
    }
    if ( !all(is.numeric(delta_draws))){
        stop("each entry in delta_draws must be numeric")
    }
    N <- dim(delta_draws)[2]
    k <- dim(delta_draws)[1]
    leng <- seq(0, N - 1)
    for (t in 1:k) {
        plot(delta_draws[t, ], xlab = '', ylab = '', type = "l", cex.main = 1.5, main = expression(paste(Trace~plot~delta)))
    }
}

#' Inefficiency Factor for Ordinal Models
#' with more than 3 outcomes
#'
#' This function calculates the inefficiency factor from the MCMC draws
#' of \eqn{(\beta, \delta)} for an ordinal model with more than 3 outcomes. The
#' inefficiency factor is calculated using the batch-means method.
#'
#' @param beta_draws      Gibbs draw of coefficients of dimension \eqn{(k x N)}.
#' @param nlags           scalar variable with default = 2.
#' @param delta_draws     Gibbs draw of cut-points.
#'
#' @details
#' Calculates the inefficiency factor of \eqn{(\beta, \delta)} using the batch-means
#' method.
#'
#' @return Returns a list with components
#' \itemize{
#' \item{\code{inefficiency_delta}: }{a vector with inefficiency factor for each \eqn{\delta}.}
#' \item{\code{inefficiency_beta}: }{a vector with inefficiency facor for each \eqn{\beta}.}
#' }
#'
#' @references Greenberg, E. (2012). “Introduction to Bayesian Econometrics.” Cambridge University
#' Press, Cambridge.
#'
#' @importFrom "pracma" "Reshape" "std"
#' @seealso pracma
#' @examples
#' set.seed(101)
#' data("data25j4")
#' x <- data25j4$x
#' y <- data25j4$y
#' p <- 0.25
#' output <- quan_regg3(y, x, mc = 50, p, 0.75)
#' beta_draws <- output$beta_draws
#' delta_draws <- output$delta_draws
#' nlags <- 2
#' inefficiency <- inefficiency_factorg3(beta_draws, nlags, delta_draws)
#'
#' # inefficiency_delta
#' #   1.433285
#' #   1.431527
#' # inefficiency_beta
#' #   1.147039
#' #   1.312946
#' #   1.423758
#'
#' @export
inefficiency_factorg3 <- function(beta_draws, nlags = 2, delta_draws) {
    if ( !all(is.numeric(beta_draws))){
        stop("each entry in beta must be numeric")
    }
    if ( !all(is.numeric(delta_draws))){
        stop("each entry in delta must be numeric")
    }
    n <- dim(beta_draws)[2]
    k <- dim(beta_draws)[1]
    inefficiency_beta <- array(0, dim = c(k, 1))
    nbatch <- floor(n / nlags)
    nuse <- nbatch * nlags
    for (i in 1:k) {
        b <- beta_draws[i, 1:nuse]
        xbatch <- Reshape(b, nlags, nbatch)
        mxbatch <- colMeans(xbatch)
        varxbatch <- sum( (t(mxbatch) - mean(b))
                          * (t(mxbatch) - mean(b))) / (nbatch - 1)
        nse <- sqrt(varxbatch / (nbatch))
        rne <- (std(b, 1) / sqrt( nuse )) / nse
        inefficiency_beta[i, 1] <- 1 / rne
    }
    n2 <- dim(delta_draws)[2]
    k2 <- dim(delta_draws)[1]
    inefficiency_delta <- array(0, dim = c(k2, 1))
    nbatch2 <- floor(n2 / nlags)
    nuse2 <- nbatch2 * nlags
    for (i in 1:k2) {
        d <- delta_draws[i, 1:nuse2]
        xbatch2 <- Reshape(d, nlags, nbatch2)
        mxbatch2 <- colMeans(xbatch2)
        varxbatch2 <- sum( (t(mxbatch2) - mean(d))
                           * (t(mxbatch2) - mean(d))) / (nbatch2 - 1)
        nse2 <- sqrt(varxbatch2 / (nbatch2))
        rne2 <- (std(d, 1) / sqrt( nuse2 )) / nse2
        inefficiency_delta[i, 1] <- 1 / rne2
    }
    result <- list("inefficiency_delta" = inefficiency_delta,
                   "inefficiency_beta" = inefficiency_beta)

    return(result)
}
