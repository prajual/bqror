#' Bayesian quantile regression for ordinal quantile model
#' with more than 3 outcomes
#'
#' This function estimates Bayesian quantile regression for ordinal quantile model with
#' more than 3 outcomes and reports the posterior mean, posterior standard deviation, and 95
#' percent posterior credible intervals of \eqn{(\beta, \delta)}.
#'
#' @usage quantregOR1(y, x, b0, B0, d0, D0, mcmc, p, tune, verbose)
#'
#' @param y                 observed ordinal outcomes, column vector of dimension \eqn{(n x 1)}.
#' @param x                 covariate matrix of dimension \eqn{(n x k)} including a column of ones with or without column names.
#' @param b0                prior mean for normal distribution to sample \eqn{\beta}, default is 0.
#' @param B0                prior variance for normal distribution to sample \eqn{\beta}.
#' @param d0                prior mean of normal distribution to sample \eqn{\delta}, default is 0.
#' @param D0                prior variance for normal distribution to sample \eqn{\delta}.
#' @param mcmc              number of MCMC iterations, post burn-in, default is 15000.
#' @param p                 quantile level or skewness parameter, p in (0,1).
#' @param tune              tuning parameter to adjust MH acceptance rate, default is 0.1.
#' @param verbose           whether to print the final output and provide additional information or not, default is TRUE.
#'
#' @details
#' Function implements the Bayesian quantile regression for
#' ordinal model with more than 3 outcomes using a combination of Gibbs sampling
#' and Metropolis-Hastings algorithm.
#'
#' Function initializes prior and then iteratively
#' samples \eqn{\beta}, \eqn{\delta} and latent variable z.
#' Burn-in is taken as \eqn{0.25*mcmc} and \eqn{nsim = burn}-\eqn{in + mcmc}.
#'
#' @return Returns a list with components:
#' \itemize{
#' \item{\code{summary}: }{summary of the MCMC draws.}
#' \item{\code{postMeanbeta}: }{vector with mean of sampled
#'  \eqn{\beta} for each covariate.}
#'  \item{\code{postMeandelta}: }{vector with mean of sampled
#'  \eqn{\delta} for each cut point.}
#'  \item{\code{postStdbeta}: }{vector with standard deviation
#'  of sampled \eqn{\beta} for each covariate.}
#'  \item{\code{postStddelta}: }{vector with standard deviation
#'  of sampled \eqn{\delta} for each cut point.}
#'  \item{\code{gamma}: }{vector of cut points including Inf and
#'  -Inf.}
#'  \item{\code{catt}}
#'  \item{\code{acceptancerate}: }{scalar to judge the acceptance
#'  rate of samples.}
#'  \item{\code{allQuantDIC}: }{results of the DIC criteria.}
#'  \item{\code{logMargLikelihood}: }{scalar value for log marginal likelihood.}
#'  \item{\code{beta}: }{matrix with all sampled values for \eqn{\beta}.}
#'  \item{\code{delta}: }{matrix with all sampled values for \eqn{\delta}.}
#' }
#'
#' @references Rahman, M. A. (2016). “Bayesian
#' Quantile Regression for Ordinal Models.”
#' Bayesian Analysis, 11(1): 1-24. DOI: 10.1214/15-BA939
#'
#' Yu, K., and Moyeed, R. A. (2001). “Bayesian Quantile Regression.” Statistics and
#' Probability Letters, 54(4): 437–447. DOI: 10.12691/ajams-6-6-4
#'
#' Casella, G., and George, E. I. (1992). “Explaining the Gibbs Sampler.”
#' The American Statistician, 46(3): 167-174. DOI: 10.1080/00031305.1992.10475878
#'
#' Geman, S., and Geman, D. (1984). “Stochastic Relaxation,
#' Gibbs Distributions, and the Bayesian Restoration of Images.”
#' IEEE Transactions an Pattern Analysis and Machine Intelligence,
#' 6(6): 721-741. DOI:10.1109/TPAMI.1984.4767596
#'
#' Chib, S., and Greenberg, E. (1995). “Understanding the Metropolis-Hastings
#' Algorithm.” The American Statistician, 49(4): 327-335. DOI: 10.2307/2684568
#'
#' Hastings, W. K. (1970). “Monte Carlo Sampling Methods Using
#' Markov Chains and Their Applications.” Biometrika, 57: 1317-1340. DOI: 10.2307/1390766
#'
#' @importFrom "stats" "sd"
#' @importFrom "stats" "rnorm" "qnorm"
#' @importFrom "pracma" "inv"
#' @importFrom "progress" "progress_bar"
#' @seealso \link[stats]{rnorm}, \link[stats]{qnorm},
#' Gibbs sampler, Metropolis-Hastings algorithm
#' @examples
#'  set.seed(101)
#'  data("data25j4")
#'  x <- data25j4$x
#'  y <- data25j4$y
#'  k <- dim(x)[2]
#'  J <- dim(as.array(unique(y)))[1]
#'  D0 <- 0.25*diag(J - 2)
#'  output <- quantregOR1(y = y,x = x, b0 = 0, B0 = 10*diag(k), d0 = 0, D0 = D0,
#'  mcmc = 40, p = 0.25, tune = 1, verbose = TRUE)
#'
#'
#'  # Number of burn-in draws: 10
#'  # Number of retained draws: 40
#'  # Summary of MCMC draws:
#'
#'
#'  #             Post Mean  Post Std   Upper Credible Lower Credible
#'  # beta_0       -2.6202   0.3588        -2.0560        -3.3243
#'  # beta_1        3.1670   0.5894         4.1713         2.1423
#'  # beta_2        4.2800   0.9141         5.7142         2.8625
#'  # delta_1       0.2188   0.4043         0.6541        -0.4384
#'  # delta_2       0.4567   0.3055         0.7518        -0.2234
#'
#'  # MH acceptance rate: 36
#'  # Log of Marginal Likelihood: -554.61
#'  # DIC: 1375.33
#'
#' @export
quantregOR1 <- function(y, x, b0 = 0, B0, d0 = 0, D0, mcmc = 15000, p, tune = 0.1, verbose = TRUE) {
    cols <- colnames(x)
    names(x) <- NULL
    names(y) <- NULL
    x <- as.matrix(x)
    y <- as.matrix(y)
    if (dim(y)[2] != 1){
        stop("input y should be a column vector")
    }
    if ( any(!all(y == floor(y)))){
        stop("each entry of y must be an integer")
    }
    if ( !all(is.numeric(x))){
        stop("each entry in x must be numeric")
    }
    if ( !all(is.numeric(B0))){
        stop("each entry in B0 must be numeric")
    }
    if ( !all(is.numeric(D0))){
        stop("each entry in D0 must be numeric")
    }
    if ( ( length(mcmc) != 1) || (length(tune) != 1 )){
        if (length(mcmc) != 1){
            stop("parameter mcmc must be scalar")
        }
        else{
            stop("parameter tune must be scalar")
        }
    }
    if (!is.numeric(mcmc)){
        stop("parameter mcmc must be a numeric")
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
    if (any(tune < 0)){
        stop("parameter tune must be greater than 0")
    }
    J <- dim(as.array(unique(y)))[1]
    if ( J <= 3 ){
        warning("The outcome variable has only 3 outcome categories.
                It is recommended to use quantregOR2 function for modelling
                it as it is more efficient and fast")
    }
    n <- dim(x)[1]
    k <- dim(x)[2]
    if ((dim(D0)[1] != (J-2)) | (dim(D0)[2] != (J-2))){
        stop("D0 is the prior variance to sample delta
             must have dimension (J-2)x(J-2)")
    }
    if ((dim(B0)[1] != (k)) | (dim(B0)[2] != (k))){
        stop("B0 is the prior variance to sample beta
             must have dimension kxk")
    }
    burn <- 0.25 * mcmc
    nsim <- burn + mcmc

    yprob <- array(0, dim = c(n, J))
    for (i in 1:n) {
        yprob[i, y[i]] <- 1
    }
    yprob <- colSums(yprob) / n
    gam <- qnorm(cumsum(yprob[1:(J - 1)]))
    deltaIn <- t(log(gam[2:(J - 1)] - gam[1:(J - 2)]))

    b0 <- array(rep(b0, k), dim = c(k, 1))
    invB0 <- inv(B0)
    invB0b0 <- invB0 %*% b0

    d0 <- array(d0, dim = c(J-2, 1))

    beta <- array (0, dim = c(k, nsim))
    delta <- array(0, dim = c((J - 2), nsim))

    ytemp <- y - 1.5
    beta[, 1] <- mldivide( (t(x) %*% (x)), (t(x) %*% ytemp))
    delta[, 1] <- deltaIn
    w <- array( (abs(rnorm(n, mean = 2, sd = 1))), dim = c (n, 1))
    z <- array( (rnorm(n, mean = 0, sd = 1)), dim = c(n, 1))

    theta <- (1 - (2 * p)) / (p * (1 - p))
    tau <- sqrt(2 / (p * (1 - p)))
    tau2 <- tau ^ 2
    lambda <- 0.5
    cri0     <- 1;
    cri1     <- 0.001;
    stepsize <- 1;
    maxiter  <- 10;
    h        <- 0.002;
    dh       <- 0.0002;
    sw       <- 20;
    minimize <- qrminfundtheorem(deltaIn, y, x,
                                 beta[, 1], cri0, cri1,
                                 stepsize, maxiter, h, dh, sw, p)

    Dhat <- -inv(minimize$H) * 3

    mhacc <- 0
    if(verbose) {
        pb <- progress_bar$new(" Simulation in Progress [:bar] :percent",
                           total = nsim, clear = FALSE, width = 100)
    }
    for (i in 2:nsim) {
        betadraw <- drawbetaOR1(z, x, w, tau2, theta, invB0, invB0b0)
        beta[, i] <- betadraw$beta

        w <- drawwOR1(z, x, beta[, i], tau2, theta, lambda)

        deltarw <- drawdeltaOR1(y, x, beta[, i], delta[, (i - 1)], d0, D0, tune, Dhat, p)
        delta[, i] <- deltarw$deltareturn
        if (i > burn) {
            mhacc <- mhacc + deltarw$accept
        }

        z <- drawlatentOR1(y, x, beta[, i], w, theta, tau2, delta[, i])
        if(verbose) {
            pb$tick()
        }
    }

    postMeanbeta <- rowMeans(beta[, (burn + 1):nsim])
    postStdbeta <- apply(beta[, (burn + 1):nsim], 1, sd)


    if(J==3) {
        postMeandelta <- mean(delta[(burn + 1):nsim])
        postStddelta <- std(delta[(burn + 1):nsim])
    }
    else {
        postMeandelta <- rowMeans(delta[, (burn + 1):nsim])
        postStddelta <- apply(delta[, (burn + 1):nsim], 1, sd)
    }

    gammacp <- array(0, dim = c(J - 1, 1))
    expdelta <- exp(postMeandelta)
    for (j in 2:(J - 1)) {
        gammacp[j] <- sum(expdelta[1:(j - 1)])
    }

    acceptrate <- (mhacc / mcmc) * 100

    xbar <- colMeans(x)
    catt <- array(0, dim = c(J))
    catt[1] <- alcdfstd( (0 - xbar %*% postMeanbeta), p)
    for (j in 2:(J - 1)) {
        catt[j] <- alcdfstd( (gammacp[j] - xbar %*% postMeanbeta), p) -
                        alcdfstd( (gammacp[(j - 1)] - xbar %*% postMeanbeta), p)
    }
    catt[J] <- 1 - alcdfstd( (gammacp[(J - 1)] - xbar %*% postMeanbeta), p)

    allQuantDIC <- devianceOR1(y, x, delta, burn, nsim,
                             postMeanbeta, postMeandelta, beta, p)

    logMargLikelihood <- logMargLikelihoodOR1(y, x, b0, B0,
                                               d0, D0, postMeanbeta,
                                               postMeandelta, beta, delta,
                                               tune, Dhat, p, verbose)

    postMeanbeta <- array(postMeanbeta, dim = c(k, 1))
    postStdbeta <- array(postStdbeta, dim = c(k, 1))
    postMeandelta <- array(postMeandelta, dim = c(J-2, 1))
    postStddelta <- array(postStddelta, dim = c(J-2, 1))

    upperCrediblebeta <- array(apply(beta[, (burn + 1):nsim], 1, quantile, c(0.975)), dim = c(k, 1))
    lowerCrediblebeta <- array(apply(beta[, (burn + 1):nsim], 1, quantile, c(0.025)), dim = c(k, 1))

    if(J==3) {
        upperCredibledelta <- quantile(delta[(burn + 1):nsim], c(0.975))
        lowerCredibledelta <- quantile(delta[(burn + 1):nsim], c(0.025))
    }
    else {
        upperCredibledelta <- array(apply(delta[, (burn + 1):nsim], 1, quantile, c(0.975)), dim = c(J-2, 1))
        lowerCredibledelta <- array(apply(delta[, (burn + 1):nsim], 1, quantile, c(0.025)), dim = c(J-2, 1))
    }

    allQuantbeta <- cbind(postMeanbeta, postStdbeta, upperCrediblebeta, lowerCrediblebeta)
    allQuantdelta <- cbind(postMeandelta, postStddelta, upperCredibledelta, lowerCredibledelta)
    summary <- rbind(allQuantbeta, allQuantdelta)
    name <- list('Post Mean', 'Post Std', 'Upper Credible', 'Lower Credible')
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
    if (verbose) {
        print(noquote(paste0("Number of burn-in draws: ", burn)))
        print(noquote(paste0("Number of retained draws: ", mcmc)))
        print(noquote("Summary of MCMC draws: "))
        cat("\n")
        print(round(summary, 4))
        cat("\n")
        print(noquote(paste0("MH acceptance rate: ", round(acceptrate, 2))))
        print(noquote(paste0('Log of Marginal Likelihood: ', round(logMargLikelihood, 2))))
        print(noquote(paste0("DIC: ", round(allQuantDIC$DIC, 2))))
    }

    result <- list("summary" = round(summary, 4),
                   "postMeanbeta" = postMeanbeta,
                   "postMeandelta" = postMeandelta,
                   "postStdbeta" = postStdbeta,
                   "postStddelta" = postStddelta,
                   "gamma" = gammacp,
                   "catt" = catt,
                   "acceptancerate" = acceptrate,
                   "allQuantDIC" = allQuantDIC,
                   "logMargLikelihood" = logMargLikelihood,
                   "beta" = beta,
                   "delta" = delta)
    return(result)
}
#' Minimize the negative of log-likelihood
#'
#' This function minimizes the negative of the log-likelihood for
#' ordinal quantile model with respect to cut-points \eqn{\delta} using the
#' fundamental theorem of calculus.
#'
#' @usage qrminfundtheorem(deltaIn, y, x, beta, cri0, cri1, stepsize, maxiter, h, dh, sw, p)
#'
#' @param deltaIn   initialization of cut-points.
#' @param y         observed ordinal outcomes, column vector of dimension \eqn{(n x 1)}.
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
#' Second derivative from first principle
#'
#' \deqn{f'(x-h)=(f(x)-f(x-h))/h}
#' \deqn{f''(x)= [{(f(x+h)-f(x))/h} - (f(x)-f(x-h))/h]/h}
#'       \deqn{= [(f(x+h)+f(x-h)-2 f(x))]/h^2}
#'
#' cross partial derivatives
#'
#' \deqn{f(x) = [f(x+dh,y)-f(x-dh,y)]/2dh}
#' \deqn{f(x,y)=[{(f(x+dh,y+dh) - f(x+dh,y-dh))/2dh} - {(f(x-dh,y+dh) -
#' f(x-dh,y-dh))/2dh}]/2dh}
#' \deqn{= 0.25* [{(f(x+dh,y+dh)-f(x+dh,y-dh))} -{(f(x-dh,y+dh)-f(x-dh,y-dh))}]/dh2}
#'
#' @return Returns a list with components
#' \itemize{
#' \item{\code{deltamin}: }{vector with cutpoints that minimize the log-likelihood function.}
#' \item{\code{negsum}: }{scalar with sum of log-likelihood values.}
#' \item{\code{logl}: }{vector with log-likelihood values.}
#' \item{\code{G}: }{gradient vector, \eqn{(n x k)} matrix with i-th row as the score
#' for the i-th unit.}
#' \item{\code{H}: }{represents Hessian matrix.}
#' }
#'
#' @importFrom "pracma" "mldivide"
#' @references Rahman, M. A. (2016). “Bayesian
#' Quantile Regression for Ordinal Models.”
#' Bayesian Analysis, 11(1): 1-24. DOI: 10.1214/15-BA939
#'
#'
#' @seealso differential calculus, functional maximization,
#' \link[pracma]{mldivide}
#' @examples
#' set.seed(101)
#' deltaIn <- c(-0.002570995,  1.044481071)
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
#' output <- qrminfundtheorem(deltaIn, y, x, beta, cri0, cri1, stepsize, maxiter, h, dh, sw, p)
#'
#' # deltamin
#' #   0.8266967 0.3635708
#' # negsum
#' #   645.4911
#' # logl
#' #    -0.7136999
#' #    -1.5340787
#' #    -1.1072447
#' #    -1.4423124
#' #    -1.3944677
#' #    -0.7941271
#' #    -1.6544072
#' #    -0.3246632
#' #    -1.8582422
#' #    -0.9220822
#' #    -2.1117739 .. soon
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
#' #   -338.21243  -41.10775
#' #   -41.10775 -106.32758
#'
#' @export
qrminfundtheorem <- function(deltaIn, y, x, beta, cri0, cri1,
                             stepsize, maxiter, h, dh, sw, p) {
    if ( !all(is.numeric(deltaIn))){
        stop("each entry in deltaIn must be numeric")
    }
    if (dim(y)[2] != 1){
        stop("input y should be a column vector")
    }
    if ( any(!all(y == floor(y)))){
        stop("each entry of y must be an integer")
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
    d <- length(deltaIn)
    storevn <- array(0, dim = c(n, d))
    storevp <- array(0, dim = c(n, d))
    checkoutput <- array(0, dim = c(n, 3 * d + 1))

    cri <- cri0
    der <- array(0, dim = c(n, d))
    dh2 <- dh ^ 2

    jj <- 0
    while ( ( cri > cri1 ) && ( jj < maxiter )) {
        jj <- jj + 1
        Quantity <- qrnegLogLikensumOR1(deltaIn, y, x, beta, p)
        vo <- -Quantity$nlogl
        deltao <- deltaIn
        for (i in 1:d) {
            deltaIn[i] <- deltaIn[i] - h
            Quantity1 <- qrnegLogLikensumOR1(deltaIn, y, x, beta, p)
            vn <- -Quantity1$nlogl
            deltaIn <- deltao

            storevn[, i] <- vn
            deltaIn[i] <- deltaIn[i] + h
            Quantity2 <- qrnegLogLikensumOR1(deltaIn, y, x, beta, p)
            vp <- -Quantity2$nlogl
            deltaIn <- deltao

            storevp[, i] <- vp
            der[, i] <- ( (0.5 * (vp - vn))) / h;
        }
        hess <- array(0, dim = c(d, d))
        i <- 1
        j <- 1
        while (j <= d) {
            while (i <= j) {
                if (i == j) {
                    deltaIn[i] <- deltaIn[i] + dh
                    Quantity3 <- qrnegLogLikensumOR1(deltaIn, y, x, beta, p)
                    vp2 <- -Quantity3$nlogl
                    deltaIn <- deltao

                    deltaIn[i] <- deltaIn[i] - dh
                    Quantity4 <- qrnegLogLikensumOR1(deltaIn, y, x, beta, p)
                    vn2 <- -Quantity4$nlogl
                    deltaIn <- deltao

                    hess[i, j] <- sum( (vp2 + vn2 - (2 * vo)) / dh2)
                }
                else{
                    f <- c(i, j)
                    deltaIn[f] <- deltaIn[f] + dh
                    Quantity5 <- qrnegLogLikensumOR1(deltaIn, y, x, beta, p)
                    vpp <- -Quantity5$nlogl
                    deltaIn <- deltao

                    deltaIn[f] <- deltaIn[f] - dh
                    Quantity6 <- qrnegLogLikensumOR1(deltaIn, y, x, beta, p)
                    vnn <- -Quantity6$nlogl
                    deltaIn <- deltao

                    deltaIn[f] <- deltaIn[f] + c(dh, -dh)
                    Quantity7 <- qrnegLogLikensumOR1(deltaIn, y, x, beta, p)
                    vpn <- -Quantity7$nlogl
                    deltaIn <- deltao

                    deltaIn[f] <- deltaIn[f] + c(-dh, dh)
                    Quantity8 <- qrnegLogLikensumOR1(deltaIn, y, x, beta, p)
                    vnp <- -Quantity8$nlogl
                    deltaIn <- deltao

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
        ddeltabhhh <- mldivide(t(der) %*% der, (colSums(der)))
        ddeltahess <- mldivide(-hess, (colSums(der)))

        ddelta <- ( ( (1 - min(1, max(0, jj - sw))) * ddeltabhhh) +
                        (min(1, max(0, jj - sw)) * ddeltahess))
        deltaIn <- deltaIn + stepsize * t(ddelta)

        if (jj == maxiter){
        }
    }
    deltamin <- deltaIn
    Quantity9 <- qrnegLogLikensumOR1(deltamin, y, x, beta, p)
    logl <- -Quantity9$nlogl
    negsum <- Quantity9$negsumlogl
    G <- der
    H <- hess
    rt <- list("deltamin" = deltamin,
               "negsum" = negsum,
               "logl" = logl,
               "G" = G,
               "H" = H)
    return(rt)
}

#' Negative log-likelihood for ordinal quantile model with more than 3 outcomes
#'
#' Function for calculating the negative log-likelihood for ordinal quantile model with
#' more than 3 outcomes.
#'
#' @usage qrnegLogLikensumOR1(deltaIn, y, x, beta, p)
#'
#' @param deltaIn   initialization of cut-points.
#' @param y         observed ordinal outcomes, column vector of dimension \eqn{(n x 1)}.
#' @param x         covariate matrix of dimension \eqn{(n x k)} including a column of ones.
#' @param beta      column vector of coefficients of dimension \eqn{(k x 1)}.
#' @param p         quantile level or skewness parameter, p in (0,1).
#'
#' @details
#' Computes the negative of the log-likelihood function for
#' ordinal quantile regression model with more than 3 outcomes.
#'
#' @return Returns a list with components
#' \itemize{
#' \item{\code{nlogl}: }{vector with likelihood values.}
#' \item{\code{negsumlogl}: }{scalar with value of negative log-likelihood.}
#' }
#'
#' @references Rahman, M. A. (2016). “Bayesian
#' Quantile Regression for Ordinal Models.”
#' Bayesian Analysis, 11(1): 1-24. DOI: 10.1214/15-BA939
#'
#' @seealso likelihood maximization
#' @examples
#' set.seed(101)
#' deltaIn <- c(-0.002570995, 1.044481071)
#' data("data25j4")
#' x <- data25j4$x
#' y <- data25j4$y
#' p <- 0.25
#' beta <- c(0.3990094, 0.8168991, 2.8034963)
#' output <- qrnegLogLikensumOR1(deltaIn, y, x, beta, p)
#'
#' # nlogl
#' #   0.7424858
#' #   1.1649645
#' #   2.1344390
#' #   0.9881085
#' #   2.7677386
#' #   0.8229129
#' #   0.8854911
#' #   0.3534490
#' #   1.8582422
#' #   0.9508680 .. soon
#'
#' # negsumlogl
#' #   663.5475
#'
#' @export
qrnegLogLikensumOR1 <- function(deltaIn, y, x, beta, p) {
    if ( !all(is.numeric(deltaIn))){
        stop("each entry in deltaIn must be numeric")
    }
    if (dim(y)[2] != 1){
        stop("input y should be a column vector")
    }
    if ( any(!all(y == floor(y)))){
        stop("each entry of y must be an integer")
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
    expdelta <- exp(deltaIn)
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
    negsumlogl <- -sum(lnpdf)
    respon <- list("nlogl" = nlogl,
                   "negsumlogl" = negsumlogl)
    return(respon)
}
#' Samples \eqn{\beta} for ordinal quantile model
#' with more than 3 outcomes
#'
#' This function samples \eqn{\beta} from its conditional
#' posterior distribution for ordinal quantile model with more than 3
#' outcomes i.e. ORI model.
#'
#' @usage drawbetaOR1(z, x, w, tau2, theta, invB0, invB0b0)
#'
#' @param z         dependent variable i.e. ordinal outcome values.
#' @param x         covariate matrix of dimension \eqn{(n x k)} including a column of ones.
#' @param w         latent weights, column vector.
#' @param tau2      2/(p(1-p)).
#' @param theta     (1-2p)/(p(1-p)).
#' @param invB0     inverse of prior covariance matrix of normal distribution.
#' @param invB0b0   prior mean pre-multiplied by invB0.
#'
#' @details
#' Function samples a vector of \eqn{\beta} from a multivariate
#' normal distribution.
#'
#' @return Returns a list with components
#' \itemize{
#' \item{\code{beta}: }{column vector of \eqn{\beta}
#' from a multivariate normal distribution.}
#' \item{\code{Btilde}: }{variance parameter for the normal
#'  distribution.}
#' \item{\code{btilde}: }{mean parameter for the
#' normal distribution.}
#' }
#'
#' @references Rahman, M. A. (2016). “Bayesian
#' Quantile Regression for Ordinal Models.”
#' Bayesian Analysis, 11(1): 1-24. DOI: 10.1214/15-BA939
#'
#' Casella, G., and George, E. I. (1992). “Explaining the Gibbs Sampler.”
#' The American Statistician, 46(3): 167-174. DOI: 10.1080/00031305.1992.10475878
#'
#' Geman, S., and Geman, D. (1984). “Stochastic Relaxation,
#' Gibbs Distributions, and the Bayesian Restoration of Images.”
#' IEEE Transactions an Pattern Analysis and Machine Intelligence,
#' 6(6): 721-741. DOI: 10.1109/TPAMI.1984.4767596
#'
#' @importFrom "MASS" "mvrnorm"
#' @importFrom "pracma" "inv"
#' @seealso Gibbs sampling, normal distribution,
#' \link[MASS]{mvrnorm}, \link[pracma]{inv}
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
#' output <- drawbetaOR1(z, x, w, tau2, theta, invB0, invB0b0)
#'
#' # output$beta
#' #   -0.2481837 0.7837995 -3.4680418
#' @export
drawbetaOR1 <- function(z, x, w, tau2, theta, invB0, invB0b0) {
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
        stop("parameter tau2 must be numeric")
    }
    if ( length(theta) != 1){
        stop("parameter theta must be scalar")
    }
    if ( !all(is.numeric(theta))){
        stop("parameter theta must be numeric")
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
    for (i in 1:n){
        meancomp[i, ] <- (x[i, ] * (z[i] - (theta * w[i]))) / (tau2 * w[i])
        varcomp[, , i] <- ( (x[i, ]) %*% (t(x[i, ]))) / (tau2 * w[i])
    }
    Btilde <- inv(invB0 + rowSums(varcomp, dims = 2))
    btilde <- Btilde %*% (invB0b0 + (colSums(meancomp)))
    L <- t(chol(Btilde))
    beta <- btilde + L %*% (mvrnorm(n = 1, mu = q, Sigma = eye))

    betaReturns <- list("beta" = beta,
                        "Btilde" = Btilde,
                        "btilde" = btilde)
    return(betaReturns)
}
#' Samples the latent weight w for ordinal quantile model
#' with more than 3 outcomes
#'
#' This function samples the latent weight w from a generalized
#' inverse-Gaussian distribution (GIG) for ordinal quantile model with more
#' than 3 outcomes.
#'
#' @usage drawwOR1(z, x, beta, tau2, theta, lambda)
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
#' @return column vector of w from a GIG distribution.
#'
#' @references Albert, J., and Chib, S. (1993). “Bayesian Analysis of Binary and Polychotomous
#' Response Data.” Journal of the American Statistical
#' Association, 88(422): 669–679. DOI: 10.1080/01621459.1993.10476321
#'
#' Casella, G., and George, E. I. (1992). “Explaining the Gibbs Sampler.”
#' The American Statistician, 46(3): 167-174. DOI: 10.1080/00031305.1992.10475878
#'
#' Geman, S., and Geman, D. (1984). “Stochastic Relaxation,
#' Gibbs Distributions, and the Bayesian Restoration of Images.”
#' IEEE Transactions an Pattern Analysis and Machine Intelligence,
#' 6(6): 721-741. DOI: 10.1109/TPAMI.1984.4767596
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
#' output <- drawwOR1(z, x, beta, tau2, theta, lambda)
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
drawwOR1 <- function(z, x, beta, tau2, theta, lambda) {
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
        stop("parameter tau2 must be numeric")
    }
    if ( length(theta) != 1){
        stop("parameter theta must be scalar")
    }
    if ( !all(is.numeric(theta))){
        stop("parameter theta must be numeric")
    }
    if ( length(lambda) != 1){
        stop("parameter lambda must be scalar")
    }
    if ( !all(is.numeric(lambda))){
        stop("parameter lambda must be numeric")
    }
    n <- dim(x)[1]
    tildeeta2 <- ( (theta ^ 2) / (tau2)) + 2
    tildelambda2 <- array(0, dim = c(n, 1))
    w <- array(0, dim = c(n, 1))
    for (i in 1:n) {
        tildelambda2[i, 1] <- ( (z[i] - (x[i, ] %*% beta)) ^ 2) / (tau2)
        w[i, 1] <- rgig(n = 1, lambda = lambda,
                        chi = tildelambda2[i, 1],
                        psi = tildeeta2)
    }
    return(w)
}
#' Samples the latent variable z for ordinal quantile model
#' with more than 3 outcomes
#'
#' This function samples the latent variable z from a truncated
#' normal distribution for ordinal quantile model with more than 3 outcomes.
#'
#' @usage drawlatentOR1(y, x, beta, w, theta, tau2, delta)
#'
#' @param y         observed ordinal outcomes, column vector of dimension \eqn{(n x 1)}.
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
#' @return column vector of values for latent variable, z.
#'
#' @references Albert, J., and Chib, S. (1993). “Bayesian Analysis of Binary and Polychotomous
#' Response Data.” Journal of the American Statistical
#' Association, 88(422): 669–679. DOI: 10.1080/01621459.1993.10476321
#'
#' Casella, G., and George, E. I. (1992). “Explaining the Gibbs Sampler.”
#' The American Statistician, 46(3): 167-174. DOI: 10.1080/00031305.1992.10475878
#'
#' Geman, S., and Geman, D. (1984). “Stochastic Relaxation,
#' Gibbs Distributions, and the Bayesian Restoration of Images.”
#' IEEE Transactions an Pattern Analysis and Machine Intelligence,
#' 6(6): 721-741. DOI: 10.1109/TPAMI.1984.4767596
#'
#' Robert, C. P. (1995). “Simulation of truncated normal variables.” Statistics and
#' Computing, 5: 121–125. DOI: 10.1007/BF00143942
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
#' output <- drawlatentOR1(y, x, beta, w, theta, tau2, delta)
#'
#' # output
#' #   0.6261896 3.129285 2.659578 8.680291
#' #   13.22584 2.545938 1.507739 2.167358
#' #   15.03059 -3.963201 9.237466 -1.813652
#' #   2.718623 -3.515609 8.352259 -0.3880043
#' #   -0.8917078 12.81702 -0.2009296 1.069133 ... soon
#'
#' @export
drawlatentOR1 <- function(y, x, beta, w, theta, tau2, delta) {
    if (dim(y)[2] != 1){
        stop("input y should be a column vector")
    }
    if ( any(!all(y == floor(y)))){
        stop("each entry of y must be an integer")
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
        stop("parameter tau2 must be numeric")
    }
    if ( length(theta) != 1){
        stop("parameter theta must be scalar")
    }
    if ( !all(is.numeric(theta))){
        stop("parameter theta must be numeric")
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
        meanp <- (x[i, ] %*% beta) + (theta * w[i])
        std <- sqrt(tau2 * w[i])
        temp <- y[i]
        a1 <- gammacp[temp]
        b1 <- gammacp[temp + 1]
        z[1, i] <- rtruncnorm(n = 1, a = a1, b = b1, mean = meanp, sd = std)
    }
    return(z)
}
#' Samples \eqn{\delta} for ordinal quantile model
#' with more than 3 outcomes
#'
#' This function samples the \eqn{\delta} using a
#' random-walk Metropolis-Hastings algorithm for ordinal
#' quantile model with more than 3 outcomes.
#'
#' @usage drawdeltaOR1(y, x, beta, delta0, d0, D0, tune, Dhat, p)
#'
#' @param y         observed ordinal outcomes, column vector of dimension \eqn{(n x 1)}.
#' @param x         covariate matrix of dimension \eqn{(n x k)} including a column of ones.
#' @param beta      Gibbs draw of coefficients of dimension \eqn{(k x 1)}.
#' @param delta0    initial value for \eqn{\delta}.
#' @param d0        prior mean of normal distribution.
#' @param D0        prior variance for normal distribution to sample \eqn{\delta}.
#' @param tune      tuning parameter to adjust MH acceptance rate.
#' @param Dhat      negative inverse Hessian from maximization of log-likelihood.
#' @param p         quantile level or skewness parameter, p in (0,1).
#'
#' @details
#' Samples the \eqn{\delta} using a random-walk Metropolis-Hastings algorithm.
#'
#' @return Returns a list with components
#' \itemize{
#'  \item{\code{deltaReturn}: }{vector with \eqn{\delta} values using MH algorithm.}
#'   \item{\code{accept}: }{indicator for acceptance of proposed value of \eqn{\delta}.}
#' }
#'
#' @references Rahman, M. A. (2016). “Bayesian
#' Quantile Regression for Ordinal Models.”
#' Bayesian Analysis, 11(1): 1-24. DOI: 10.1214/15-BA939
#'
#' Chib, S., and Greenberg, E. (1995). “Understanding the Metropolis-Hastings
#' Algorithm.” The American Statistician, 49(4): 327-335. DOI: 10.2307/2684568
#'
#' Hastings, W. K. (1970). “Monte Carlo Sampling Methods Using
#' Markov Chains and Their Applications.” Biometrika, 57: 1317-1340. DOI: 10.2307/1390766
#'
#' Jeliazkov, I., Graves, J., and Kutzbach, M. (2008). “Fitting and Comparison of Models
#' for Multivariate Ordinal Outcomes.” Advances in Econometrics: Bayesian Econometrics,
#' 23: 115–156. DOI: 10.1016/S0731-9053(08)23004-5
#'
#' Jeliazkov, I., and Rahman, M. A. (2012). “Binary and Ordinal Data Analysis
#' in Economics: Modeling and Estimation” in Mathematical Modeling with Multidisciplinary
#' Applications, edited by X.S. Yang, 123-150. John Wiley & Sons Inc, Hoboken, New Jersey. DOI: 10.1002/9781118462706.ch6
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
#' output <- drawdeltaOR1(y, x, beta, delta0, d0, D0, tune, Dhat, p)
#'
#' # deltareturn
#' #   -0.9025802 -2.229514
#' # accept
#' #   1
#'
#' @export
drawdeltaOR1 <- function(y, x, beta, delta0, d0, D0, tune, Dhat, p){
    if (dim(y)[2] != 1){
        stop("input y should be a column vector")
    }
    if ( any(!all(y == floor(y)))){
        stop("each entry of y must be an integer")
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
    if ( !all(is.numeric(Dhat))){
        stop("each entry in Dhat must be numeric")
    }
    if (any(p < 0 | p > 1)){
        stop("parameter p must be between 0 to 1")
    }
    J <- dim(as.array(unique(y)))[1]
    k <- (J - 2)
    L <- t(chol(Dhat))
    delta1 <- delta0 + tune * t(L %*%  (rnorm(n = k, mean = 0, sd = 1)))
    num <-  qrnegLogLikensumOR1(delta1, y, x, beta, p)
    den <-  qrnegLogLikensumOR1(delta0, y, x, beta, p)
    pnum <- -num$negsumlogl +
                   mvnpdf(x = matrix(delta1),
                          mean  = matrix(d0),
                          varcovM = D0,
                          Log = TRUE)
    pden <- -den$negsumlogl +
                   mvnpdf(x = matrix(delta0),
                          mean  = matrix(d0),
                          varcovM = D0,
                          Log = TRUE)
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
#' Deviance Information Criteria for ordinal quantile model
#' with more than 3 outcomes
#'
#' Function for computing the Deviance information criteria for ordinal quantile
#' model with more than 3 outcomes.
#'
#' @usage devianceOR1(y, x, deltastore, burn, nsim, postMeanbeta, postMeandelta, beta, p)
#'
#' @param y                observed ordinal outcomes, column vector of dimension \eqn{(n x 1)}.
#' @param x                covariate matrix of dimension \eqn{(n x k)} including a column of ones.
#' @param postMeanbeta     mean value of \eqn{\beta} obtained from MCMC draws.
#' @param postMeandelta    mean value of \eqn{\delta} obtained from MCMC draws.
#' @param beta             MCMC draw of coefficients, dimension is \eqn{(k x nsim)}.
#' @param deltastore       MCMC draws of \eqn{\delta}.
#' @param p                quantile level or skewness parameter, p in (0,1).
#' @param burn             number of discarded MCMC iterations.
#' @param nsim             total number of samples, including the burn-in.
#'
#' @details
#' Deviance is -2*(log likelihood) and has an important role in
#' statistical model comparison because of its relation with Kullback-Leibler
#' information criteria.
#'
#' @return Returns a list with components
#' \deqn{DIC = 2*avgdDeviance - devpostmean}
#' \deqn{pd = avgdDeviance - devpostmean}
#' \deqn{devpostmean = -2*(logLikelihood)}.
#'
#' @references Rahman, M. A. (2016). “Bayesian
#' Quantile Regression for Ordinal Models.”
#' Bayesian Analysis, 11(1): 1-24. DOI: 10.1214/15-BA939
#'
#' Spiegelhalter, D. J., Best, N. G., Carlin, B. P. and Linde, A. (2002).
#' “Bayesian Measures of Model Complexity and Fit.” Journal of the
#' Royal Statistical Society B, Part 4: 583-639. DOI: 10.1111/1467-9868.00353
#'
#' Gelman, A., Carlin, J. B., Stern, H. S., and Rubin, D. B.
#' “Bayesian Data Analysis.” 2nd Edition, Chapman and Hall. DOI: 10.1002/sim.1856
#'
#' @seealso  decision criteria
#' @examples
#' set.seed(101)
#' data("data25j4")
#' x <- data25j4$x
#' y <- data25j4$y
#' k <- dim(x)[2]
#' J <- dim(as.array(unique(y)))[1]
#' D0 <- 0.25*diag(J - 2)
#' output <- quantregOR1(y = y, x = x, b0 = 0, B0 = 10*diag(k), d0 = 0, D0 = D0,
#' mcmc = 40, p = 0.25, tune = 1, verbose = FALSE)
#' mcmc <- 40
#' deltastore <- output$delta
#' burn <- 0.25*mcmc
#' nsim <- burn + mcmc
#' postMeanbeta <- output$postMeanbeta
#' postMeandelta <- output$postMeandelta
#' beta <- output$beta
#' deviance <- devianceOR1(y, x, deltastore, burn, nsim,
#' postMeanbeta, postMeandelta, beta, p = 0.25)
#'
#' # DIC
#' #   1375.329
#' # pd
#' #   139.1751
#' # devpostmean
#' #   1096.979
#'
#' @export
devianceOR1 <- function(y, x, deltastore, burn,
                       nsim, postMeanbeta, postMeandelta, beta, p){
    cols <- colnames(x)
    names(x) <- NULL
    names(y) <- NULL
    x <- as.matrix(x)
    y <- as.matrix(y)
    if (dim(y)[2] != 1){
        stop("input y should be a column vector")
    }
    if ( any(!all(y == floor(y)))){
        stop("each entry of y must be an integer")
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
    if ( length(nsim) != 1){
        stop("parameter nsim must be scalar")
    }
    if ( !all(is.numeric(postMeanbeta))){
        stop("each entry in postMeanbeta must be numeric")
    }
    if ( !all(is.numeric(postMeandelta))){
        stop("each entry in postMeandelta must be numeric")
    }
    if ( !all(is.numeric(beta))){
        stop("each entry in beta must be numeric")
    }
    if (any(p < 0 | p > 1)){
        stop("parameter p must be between 0 to 1")
    }
    delta <- deltastore
    devpostmean <- array(0, dim = c(1))
    DIC <- array(0, dim = c(1))
    pd <- array(0, dim = c(1))
    ans <- qrnegLogLikensumOR1(postMeandelta, y, x,
                            postMeanbeta, p)
    devpostmean <- 2 * ans$negsumlogl
    nsim <- dim(beta[, (burn + 1):nsim])[1]
    Deviance <- array(0, dim = c(nsim, 1))
    for (i in 1:nsim) {
        temp <- qrnegLogLikensumOR1(delta[, (burn + i)],
                                 y, x, beta[, (burn + i)], p)
        Deviance[i, 1] <- 2 * temp$negsumlogl
    }
    avgdDeviance <- mean(Deviance)
    DIC <- 2 * avgdDeviance - devpostmean
    pd <- avgdDeviance - devpostmean
    result <- list("DIC" = DIC,
                   "pd" = pd,
                   "devpostmean" = devpostmean)
    return(result)
}
#' CDF of standard asymmetric Laplace distribution
#'
#' This function computes the CDF of standard asymmetric
#' Laplace distribution i.e. AL\eqn{(0, 1 ,p)}.
#'
#' @usage alcdfstd(x, p)
#'
#' @param x     scalar value.
#' @param p     quantile level or skewness parameter, p in (0,1).
#'
#' @details
#' Computes the CDF of a standard asymmetric Laplace distribution.
#' \deqn{CDF(x) = F(x) = P(X \le x)} where X is a
#' random variable that follows AL\eqn{(0, 1 ,p)}.
#'
#' @return Returns the cumulative probability value from the CDF of an asymmetric
#' Laplace distribution.
#'
#' @references Rahman, M. A. (2016). “Bayesian
#' Quantile Regression for Ordinal Models.”
#' Bayesian Analysis, 11(1): 1-24. DOI: 10.1214/15-BA939
#'
#' Yu, K., and Zhang, J. (2005). “A Three-Parameter Asymmetric
#' Laplace Distribution.” Communications in Statistics - Theory and Methods, 34(9-10), 1867-1879. DOI: 10.1080/03610920500199018
#'
#' @seealso asymmetric Laplace distribution
#' @examples
#' set.seed(101)
#' x <-  -0.5428573
#' p <- 0.25
#' output <- alcdfstd(x, p)
#'
#' # output
#' #   0.1663873
#'
#' @export
alcdfstd <- function(x, p) {
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
#' Asymmetric Laplace distribution
#'
#' This function computes the cumulative distribution function (CDF) of
#' an asymmetric Laplace distribution.
#'
#' @usage alcdf(x, mu, sigma, p)
#'
#' @param x     scalar value.
#' @param mu    location parameter of ALD.
#' @param sigma scale parameter of ALD.
#' @param p     quantile or skewness parameter, p in (0,1).
#'
#' @details
#' Computes the cumulative distribution function of
#' the asymmetric Laplace distribution.
#' \deqn{CDF(x) = F(x) = P(X \le x)} where X is a
#' random variable
#'
#' @return Returns a scalar with cumulative probability value at
#' point “x”.
#'
#' @references Rahman, M. A. (2016). “Bayesian
#' Quantile Regression for Ordinal Models.”
#' Bayesian Analysis, 11(1): 1-24. DOI: 10.1214/15-BA939
#'
#' Yu, K., and Zhang, J. (2005). “A Three-Parameter Asymmetric
#' Laplace Distribution.” Communications in Statistics - Theory and Methods, 34(9-10), 1867-1879. DOI: 10.1080/03610920500199018
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

#' Inefficiency factor for ordinal quantile model
#' with more than 3 outcomes
#'
#' This function calculates the inefficiency factor from the MCMC draws
#' of \eqn{(\beta, \delta)} for ordinal quantile model with more than 3 outcomes. The
#' inefficiency factor is calculated using the batch-means method.
#'
#' @usage infactorOR1(x, beta, delta, autocorrelationCutoff)
#'
#' @param x                         covariate matrix of dimension \eqn{(n x k)} including a column of ones with or without column names.
#' @param beta                      Gibbs draw of coefficients of dimension \eqn{(k x nsim)}.
#' @param delta                     Gibbs draw of cut-points.
#' @param autocorrelationCutoff     cut-off to identify the number of lags, default is 0.05.
#'
#' @details
#' Calculates the inefficiency factor of \eqn{(\beta, \delta)} using the batch-means
#' method.
#'
#' @return Returns a list with components
#' \itemize{
#' \item{\code{inefficiencyDelta}: }{vector with inefficiency factor for each \eqn{\delta}.}
#' \item{\code{inefficiencyBeta}: }{vector with inefficiency factor for each \eqn{\beta}.}
#' }
#'
#' @references Greenberg, E. (2012). “Introduction to Bayesian Econometrics.” Cambridge University
#' Press, Cambridge. DOI: 10.1017/CBO9780511808920
#'
#' @importFrom "pracma" "Reshape" "std"
#' @importFrom "stats" "acf"
#' @seealso pracma, \link[stats]{acf}
#' @examples
#' set.seed(101)
#' data("data25j4")
#' x <- data25j4$x
#' y <- data25j4$y
#' k <- dim(x)[2]
#' J <- dim(as.array(unique(y)))[1]
#' D0 <- 0.25*diag(J - 2)
#' output <- quantregOR1(y = y, x = x, b0 = 0, B0 = 10*diag(k), d0 = 0, D0 = D0,
#' mcmc = 40, p = 0.25, tune = 1, verbose = FALSE)
#' beta <- output$beta
#' delta <- output$delta
#' inefficiency <- infactorOR1(x, beta, delta, 0.5)
#'
#' # Summary of Inefficiency Factor:
#'
#' #             Inefficiency
#' # beta_0        1.1008
#' # beta_1        3.0024
#' # beta_2        2.8543
#' # delta_1       3.6507
#' # delta_2       3.1784
#'
#' @export
infactorOR1 <- function(x, beta, delta, autocorrelationCutoff = 0.05) {
    cols <- colnames(x)
    names(x) <- NULL
    x <- as.matrix(x)
    if ( !all(is.numeric(beta))){
        stop("each entry in beta must be numeric")
    }
    if ( !all(is.numeric(delta))){
        stop("each entry in delta must be numeric")
    }
    n <- dim(beta)[2]
    k <- dim(beta)[1]
    inefficiencyBeta <- array(0, dim = c(k, 1))
    for (i in 1:k) {
        autocorrelation <- acf(beta[i,], plot = FALSE)
        nlags <- min(which(autocorrelation$acf <= autocorrelationCutoff))
        nbatch <- floor(n / nlags)
        nuse <- nbatch * nlags
        b <- beta[i, 1:nuse]
        xbatch <- Reshape(b, nlags, nbatch)
        mxbatch <- colMeans(xbatch)
        varxbatch <- sum( (t(mxbatch) - mean(b))
                          * (t(mxbatch) - mean(b))) / (nbatch - 1)
        nse <- sqrt(varxbatch / nbatch)
        rne <- (std(b, 1) / sqrt( nuse )) / nse
        inefficiencyBeta[i, 1] <- 1 / rne
    }
    k2 <- dim(delta)[1]
    inefficiencyDelta <- array(0, dim = c(k2, 1))
    for (i in 1:k2) {
        autocorrelation <- acf(delta[i,], plot = FALSE)
        nlags <- min(which(autocorrelation$acf <= autocorrelationCutoff))
        nbatch2 <- floor(n / nlags)
        nuse2 <- nbatch2 * nlags
        d <- delta[i, 1:nuse2]
        xbatch2 <- Reshape(d, nlags, nbatch2)
        mxbatch2 <- colMeans(xbatch2)
        varxbatch2 <- sum( (t(mxbatch2) - mean(d))
                           * (t(mxbatch2) - mean(d))) / (nbatch2 - 1)
        nse2 <- sqrt(varxbatch2 / nbatch2)
        rne2 <- (std(d, 1) / sqrt( nuse2 )) / nse2
        inefficiencyDelta[i, 1] <- 1 / rne2
    }

    inefficiencyRes <- rbind(inefficiencyBeta, inefficiencyDelta)
    name <- list('Inefficiency')
    dimnames(inefficiencyRes)[[2]] <- name
    dimnames(inefficiencyRes)[[1]] <- letters[1:(k+k2)]
    j <- 1
    if (is.null(cols)) {
        rownames(inefficiencyRes)[j] <- c('Intercept')
        for (i in paste0("beta_",1:k-1)) {
            rownames(inefficiencyRes)[j] = i
            j = j + 1
        }
    }
    else {
        for (i in cols) {
            rownames(inefficiencyRes)[j] = i
            j = j + 1
        }
    }
    for (i in paste0("delta_",1:(k2))) {
        rownames(inefficiencyRes)[j] = i
        j = j + 1
    }

    print(noquote('Summary of Inefficiency Factor: '))
    cat("\n")
    print(round(inefficiencyRes, 4))

    result <- list("inefficiencyDelta" = inefficiencyDelta,
                   "inefficiencyBeta" = inefficiencyBeta)

    return(result)
}
#' Covariate effect for Bayesian quantile regression for ordinal quantile model
#' with more than 3 outcomes
#'
#' This function computes the average covariate effect for different
#' outcomes of the ORI model at the specified quantiles. The covariate
#' effects are calculated marginally of the parameters and the remaining covariates.
#'
#' @usage covEffectOR1(model, y, x, modX, p)
#'
#' @param model     outcome of the ORI (quantregOR1) model.
#' @param y         observed ordinal outcomes, column vector of dimension \eqn{(n x 1)}.
#' @param x         covariate matrix of dimension \eqn{(n x k)} including a column of ones with or without column names.
#'                  If the covariate of interest is continuous, then the column for the covariate of interest remains unchanged.
#'                  If it is an indicator variable then replace the column for the covariate of interest with a
#'                  column of zeros.
#' @param modX      matrix x with suitable modification to an independent variable including a column of ones with
#'                  or without column names. If the covariate of interest is continuous, then add the incremental change
#'                  to each observation in the column for the covariate of interest. If the covariate is an indicator variable,
#'                  then replace the column for the covariate of interest with a column of ones.
#' @param p         quantile level or skewness parameter, p in (0,1).
#'
#' @details
#' This function computes the average covariate effect for different
#' outcomes of the ORI model at the specified quantiles. The covariate
#' effects are calculated marginally of the parameters and the remaining covariates. The computation of covariate effects utilizes
#' the MCMC outputs from estimation.
#'
#' @return Returns a list with components:
#' \itemize{
#' \item{\code{avgDiffProb}: }{vector with change in predicted
#' probabilities for each outcome category.}
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
#' Jeliazkov, I. and Rahman, M. A. (2012). “Binary and Ordinal Data Analysis
#' in Economics: Modeling and Estimation” in Mathematical Modeling with Multidisciplinary
#' Applications, edited by X.S. Yang, 123-150. John Wiley & Sons Inc, Hoboken, New Jersey. DOI: 10.1002/9781118462706.ch6
#'
#' @importFrom "stats" "sd"
#' @examples
#' set.seed(101)
#' data("data25j4")
#' x <- data25j4$x
#' y <- data25j4$y
#' k <- dim(x)[2]
#' J <- dim(as.array(unique(y)))[1]
#' D0 <- 0.25*diag(J - 2)
#' output <- quantregOR1(y = y, x = x, b0 = 0, B0 = 10*diag(k), d0 = 0, D0 = D0,
#' mcmc = 30, p = 0.25, tune = 1, verbose = FALSE)
#' modX <- x
#' modX[,3] <- modX[,3] + 0.02
#' res <- covEffectOR1(output, y, x, modX, p = 0.25)
#'
#' # Summary of Covariate Effect:
#'
#' #               Covariate Effect
#' # Category_1          -0.0076
#' # Category_2          -0.0014
#' # Category_3          -0.0010
#' # Category_4           0.0100
#'
#' @export
covEffectOR1 <- function(model, y, x, modX, p) {
    cols <- colnames(x)
    cols1 <- colnames(modX)
    names(modX) <- NULL
    names(y) <- NULL
    names(x) <- NULL
    x <- as.matrix(x)
    modX <- as.matrix(modX)
    y <- as.matrix(y)
    J <- dim(as.array(unique(y)))[1]
    if ( J <= 3 ){
        stop("This function is only available for models with more than 3 outcome
                variables.")
    }
    if (dim(y)[2] != 1){
        stop("input y should be a column vector")
    }
    if ( any(!all(y == floor(y)))){
        stop("each entry of y must be an integer")
    }
    if ( !all(is.numeric(x))){
        stop("each entry in x must be numeric")
    }
    if ( !all(is.numeric(modX))){
        stop("each entry in modX must be numeric")
    }
    if ( length(p) != 1){
        stop("parameter p must be scalar")
    }
    if (any(p < 0 | p > 1)){
        stop("parameter p must be between 0 to 1")
    }
    N <- dim(model$beta)[2]
    m <- (N)/(1.25)
    burn <- 0.25 * m
    n <- dim(x)[1]
    k <- dim(x)[2]
    betaBurnt <- model$beta[, (burn + 1):N]
    deltaBurnt <- model$delta[, (burn + 1):N]
    expdeltaBurnt <- exp(deltaBurnt)
    gammacpCov <- array(0, dim = c(J-1, m))
    for (j in 2:(J-1)) {
        gammacpCov[j, ] <- sum(expdeltaBurnt[1:(j - 1), 1])
    }
    mu <- 0
    sigma <- 1
    oldProb <- array(0, dim = c(n, m, J))
    newProb <- array(0, dim = c(n, m, J))
    oldComp <- array(0, dim = c(n, m, (J-1)))
    newComp <- array(0, dim = c(n, m, (J-1)))
    for (j in 1:(J-1)) {
        for (b in 1:m) {
            for (i in 1:n) {
                oldComp[i, b, j] <- alcdf((gammacpCov[j, b] - (x[i, ] %*% betaBurnt[, b])), mu, sigma, p)
                newComp[i, b, j] <- alcdf((gammacpCov[j, b] - (modX[i, ] %*% betaBurnt[, b])), mu, sigma, p)
            }
            if (j == 1) {
                oldProb[, b, j] <- oldComp[, b, j]
                newProb[, b, j] <- newComp[, b, j]
            }
            else {
                oldProb[, b, j] <- oldComp[, b, j] - oldComp[, b, (j-1)]
                newProb[, b, j] <- newComp[, b, j] - newComp[, b, (j-1)]
            }
        }
    }
    oldProb[, , J] = 1 - oldComp[, , (J-1)]
    newProb[, , J] = 1 - newComp[, , (J-1)]
    diffProb <- newProb - oldProb
    avgDiffProb <- array((colMeans(diffProb, dims = 2)), dim = c(J, 1))
    name <- list('Covariate Effect')
    dimnames(avgDiffProb)[[2]] <- name
    dimnames(avgDiffProb)[[1]] <- letters[1:(J)]
    ordOutput <- as.array(unique(y))
    j <- 1
    for (i in paste0("Category_",1:J)) {
        rownames(avgDiffProb)[j] = i
        j = j + 1
    }
    print(noquote('Summary of Covariate Effect: '))
    cat("\n")
    print(round(avgDiffProb, 4))

    result <- list("avgDiffProb" = avgDiffProb)

    return(result)
}
#' Marginal likelihood for ordinal quantile model
#' with more than 3 outcomes
#'
#' This function computes the logarithm of marginal likelihood for ordinal
#' quantile model with more than 3 outcomes using MCMC output from the
#' complete and reduced runs.
#'
#' @usage logMargLikelihoodOR1(y, x, b0, B0, d0, D0, postMeanbeta,
#' postMeandelta, beta, delta, tune, Dhat, p, verbose)
#'
#' @param y                 observed ordinal outcomes, column vector of dimension \eqn{(n x 1)}.
#' @param x                 covariate matrix of dimension \eqn{(n x k)} including a column of ones with or without column names.
#' @param b0                prior mean for normal distribution to sample \eqn{\beta}.
#' @param B0                prior variance for normal distribution to sample \eqn{\beta}
#' @param d0                prior mean for normal distribution to sample \eqn{\delta}.
#' @param D0                prior variance for normal distribution to sample \eqn{\delta}.
#' @param postMeanbeta      a vector with mean of sampled \eqn{\beta} for each covariate.
#' @param postMeandelta     a vector with mean of sampled \eqn{\delta} for each cut-point.
#' @param beta              a storage matrix with all sampled values for \eqn{\beta}.
#' @param delta             a storage matrix with all sampled values for \eqn{\delta}.
#' @param tune              tuning parameter to adjust MH acceptance rate.
#' @param Dhat              negative inverse Hessian from maximization of log-likelihood.
#' @param p                 quantile level or skewness parameter, p in (0,1).
#' @param verbose           whether to print the final output and provide additional information or not, default is TRUE.
#'
#' @details
#' This function computes the logarithm of marginal likelihood for ordinal quantile model with more than
#' 3 outcomes using the MCMC outputs.
#'
#' @return Returns a scalar for logarithm of marginal likelihood
#'
#' @references Rahman, M. A. (2016). “Bayesian
#' Quantile Regression for Ordinal Models.”
#' Bayesian Analysis, 11(1): 1-24. DOI: 10.1214/15-BA939
#'
#' Chib, S., and Greenberg, E. (1995). “Understanding the Metropolis-Hastings
#' Algorithm.” The American Statistician, 49(4): 327-335. DOI: 10.2307/2684568
#'
#' Chib, S. (1995). “Marginal likelihood from the Gibbs output.” Journal of the American
#' Statistical Association, 90(432):1313–1321, 1995. DOI: 10.1080/01621459.1995.10476635
#'
#' Chib, S., and Jeliazkov, I. (2001). “Marginal likelihood from the Metropolis-Hastings output.” Journal of the
#' American Statistical Association, 96(453):270–281, 2001. DOI: 10.1198/016214501750332848
#'
#' Greenberg, E. (2012). “Introduction to Bayesian Econometrics.” Cambridge University
#' Press, Cambridge. DOI: 10.1017/CBO9780511808920
#'
#' @importFrom "stats" "sd" "dnorm"
#' @importFrom "pracma" "inv"
#' @importFrom "NPflow" "mvnpdf"
#' @importFrom "progress" "progress_bar"
#'
#' @seealso \link[NPflow]{mvnpdf}, \link[stats]{dnorm},
#' Gibbs sampling, Metropolis-Hastings algorithm
#' @examples
#' set.seed(101)
#' data("data25j4")
#' x <- data25j4$x
#' y <- data25j4$y
#' k <- dim(x)[2]
#' J <- dim(as.array(unique(y)))[1]
#' D0 <- 0.25*diag(J - 2)
#' output <- quantregOR1(y = y, x = x, b0 = 0, B0 = 10*diag(k), d0 = 0, D0 = D0,
#' mcmc = 40, p = 0.25, tune = 1, verbose = FALSE)
#' # output$logMargLikelihood
#' #   -554.61
#'
#' @export
logMargLikelihoodOR1 <- function(y, x, b0, B0, d0, D0, postMeanbeta, postMeandelta, beta, delta, tune, Dhat, p, verbose) {
    cols <- colnames(x)
    names(x) <- NULL
    names(y) <- NULL
    x <- as.matrix(x)
    y <- as.matrix(y)
    if ( dim(y)[2] != 1){
        stop("input y should be a column vector")
    }
    if ( any(!all(y == floor(y)))){
        stop("each entry of y must be an integer")
    }
    if ( !all(is.numeric(x))){
        stop("each entry in x must be numeric")
    }
    if ( !all(is.numeric(B0))){
        stop("each entry in B0 must be numeric")
    }
    if ( !all(is.numeric(D0))){
        stop("each entry in D0 must be numeric")
    }
    if ( length(p) != 1){
        stop("parameter p must be scalar")
    }
    if ( any(p < 0 | p > 1)){
        stop("parameter p must be between 0 to 1")
    }
    if ( !all(is.numeric(b0))){
        stop("each entry in b0 must be numeric")
    }
    if ( !all(is.numeric(d0))){
        stop("each entry in d0 must be numeric")
    }
    if ( !all(is.numeric(Dhat))){
        stop("each entry in Dhat must be numeric")
    }
    if (any(tune < 0)){
        stop("parameter tune must be greater than 0")
    }
    J <- dim(as.array(unique(y)))[1]
    n <- dim(x)[1]
    k <- dim(x)[2]
    if ((dim(D0)[1] != (J-2)) | (dim(D0)[2] != (J-2))){
        stop("D0 is the prior variance to sample delta
             must have dimension (J-2)x(J-2)")
    }
    if ((dim(B0)[1] != (k)) | (dim(B0)[2] != (k))){
        stop("B0 is the prior variance to sample beta
             must have dimension kxk")
    }
    nsim <- dim(beta)[2]
    burn <- (0.25 * nsim) / (1.25)

    lambda <- 0.5
    theta <- (1 - 2 * p) / (p * (1 - p))
    tau <- sqrt(2 / (p * (1 - p)))
    tau2 <- tau^2
    invB0 <- inv(B0)
    invB0b0 <- invB0 %*% b0
    w <- array( (abs(rnorm(n, mean = 2, sd = 1))), dim = c (n, 1))
    z <- array( (rnorm(n, mean = 0, sd = 1)), dim = c(n, 1))

    j <- 1
    postOrddeltanum <- array(0, dim = c((nsim-burn), 1))
    postOrddeltaden <- array(0, dim = c((nsim-burn), 1))
    postOrdbetaStore <- array(0, dim = c((nsim-burn), 1))

    betaStoreRedrun <- array (0, dim = c(k, nsim))
    btildeStoreRedrun <- array(0, dim = c(k, nsim))
    BtildeStoreRedrun <- array(0, dim = c(k, k, nsim))
    deltaStoreRedrun <- array(0, dim = c((J - 2), nsim))

    if(verbose) {
        pb <- progress_bar$new(" Reduced Run in Progress [:bar] :percent",
                           total = nsim, clear = FALSE, width = 100)
    }
    for (i in 1:nsim) {
        betadrawRedrun <- drawbetaOR1(z, x, w, tau2, theta, invB0, invB0b0)

        betaStoreRedrun[, i] <- betadrawRedrun$beta
        btildeStoreRedrun[, i] <- betadrawRedrun$btilde
        BtildeStoreRedrun[, , i] <- betadrawRedrun$Btilde

        w <- drawwOR1(z, x, betaStoreRedrun[, i], tau2, theta, lambda)

        z <- drawlatentOR1(y, x, betaStoreRedrun[, i], w, theta, tau2, postMeandelta)

        deltarw <- drawdeltaOR1(y, x, betaStoreRedrun[, i], postMeandelta, d0, D0, tune, Dhat, p)
        deltaStoreRedrun[, i] <- deltarw$deltareturn
        if(verbose) {
            pb$tick()
        }
    }
    if(verbose) {
        pb <- progress_bar$new(" Calculating Marginal Likelihood [:bar] :percent",
                           total = (nsim-burn), clear = FALSE, width = 100)
    }

    for (i in (burn+1):(nsim)) {
        E1alphaMH_logLikeNum <- qrnegLogLikensumOR1(postMeandelta, y, x, beta[, i], p)
        E1alphaMH_logLikeDen <- qrnegLogLikensumOR1(delta[, i], y, x, beta[, i], p)

        E1alphaMH_logNum <- -E1alphaMH_logLikeNum$negsumlogl +
                       mvnpdf(x = matrix(postMeandelta),
                              mean  = matrix(d0),
                              varcovM = D0,
                              Log = TRUE)
        E1alphaMH_logDen <- -E1alphaMH_logLikeDen$negsumlogl +
                       mvnpdf(x = matrix(delta[, i]),
                              mean  = matrix(d0),
                              varcovM = D0,
                              Log = TRUE)
        E1alphaMH <- min(1, exp((E1alphaMH_logNum - E1alphaMH_logDen)))
        qpdf <- mvnpdf(x = matrix(postMeandelta), mean  = matrix(delta[, i]), varcovM = (tune^2)*Dhat, Log = FALSE)
        postOrddeltanum[j,] <- E1alphaMH*qpdf


        E2alphaMH_logLikeNum <- qrnegLogLikensumOR1(deltaStoreRedrun[, i], y, x, betaStoreRedrun[, i], p)
        E2alphaMH_logLikeDen <- qrnegLogLikensumOR1(postMeandelta, y, x, betaStoreRedrun[, i], p)

        E2alphaMH_logNum <- -E2alphaMH_logLikeNum$negsumlogl +
            mvnpdf(x = matrix(deltaStoreRedrun[, i]),
                   mean  = matrix(d0),
                   varcovM = D0,
                   Log = TRUE)
        E2alphaMH_logDen <- -E2alphaMH_logLikeDen$negsumlogl +
            mvnpdf(x = matrix(postMeandelta),
                   mean  = matrix(d0),
                   varcovM = D0,
                   Log = TRUE)
        postOrddeltaden[j,] <- min(1, exp((E2alphaMH_logNum - E2alphaMH_logDen)))
        postOrdbetaStore[j,] <- mvnpdf(x = matrix(postMeanbeta), mean = btildeStoreRedrun[, i], varcovM = BtildeStoreRedrun[, , i], Log = FALSE)

        j <- j + 1
        if(verbose) {
            pb$tick()
        }
    }
    postOrddelta <- mean(postOrddeltanum)/mean(postOrddeltaden)
    postOrdbeta <- mean(postOrdbetaStore)

    priorContbeta <- mvnpdf(x = matrix(postMeanbeta), mean = b0, varcovM = B0, Log = FALSE)
    priorContdelta <- mvnpdf(x = matrix(postMeandelta), mean = d0, varcovM = D0, Log = FALSE)

    logLikeCont <- -1* ((qrnegLogLikensumOR1(postMeandelta, y, x, postMeanbeta, p))$negsumlogl)
    logPriorCont <- log(priorContbeta*priorContdelta)
    logPosteriorCont <- log(postOrdbeta*postOrddelta)

    logMargLikelihood <- logLikeCont + logPriorCont - logPosteriorCont
    return(logMargLikelihood)
}
