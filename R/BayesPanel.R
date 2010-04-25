################################################################################
## This code implements the estimation of random effects panel model with the
## distributions of both the error terms and the random effects to be Gaussian.
##
## The estimation is realized through Markov chain Monte Carlo(MCMC) techniques.
## For detailed information on the methods, please refer to Wakefiled et al.
## (1994) and Chib and Carlin (1999).
##
## Codes written by Chunhua Wu and Siddhartha Chib, 04/25/2010
## For comments, please email to chunhuawu@wustl.edu
################################################################################


"BayesPanel" <-
  function(formula, data, subset, na.action, index,
           linkage = c("normal", "t"), hetero = FALSE, arma = c(0, 0),
           prior = list(beta0 = NULL, B0 = NULL, nuF = NULL, rho0 = NULL, R0 = NULL,
           nuG = NULL, nu00 = NULL, delta00 = NULL,
           nu0 = NULL,delta0 = NULL, phi0 = NULL, P0 = NULL),
           control = list(seed = 314159265, burnin = 1000, iter = 10000,
             iter2 = 10000, thin = 1, verbose = 0), ...  ){
    ## match input options
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0)
    mf <- mf[c(1, m)]
    f <- Formula(formula)
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf$formula <- f

    ## check index, sort data frame  and create ind variable
    check.index <- CheckIndex(data, index)
    data <- check.index$data
    ind <- check.index$ind
    indnames <- check.index$indnames

    mf$data <- data
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")

    y <- model.response(mf, "numeric")
    X <- model.matrix(f, data=mf, rhs=1) ##construction of X matrix
    W <- model.matrix(f, data=mf, rhs=2) ##construction of W matrix
    Xnames <- dimnames(X)[[2]]
    Wnames <- dimnames(W)[[2]]


    
    ## create dimensions
    nn <- nrow(data)
    n <- length(ind)
    k <- ncol(X)
    q <- ncol(W)

    ## check model specifications
    linkage <- match.arg(linkage)
    arma <- CheckARMA(arma)
#    ar <- arma[1]
 #   ma <- arma[2]
    ar <- 0
    ma <- 0
    
    ## check priors
    beta0 <- prior$beta0
    B0 <- prior$B0
    nuF <- prior$nuF
    rho0 <- prior$rho0
    R0 <- prior$R0
    nuG <- prior$nuG
    nu00 <- prior$nu00
    delta00 <- prior$delta00
    nu0 <- prior$nu0
    delta0 <- prior$delta0
    phi0 <- prior$phi0
    P0 <- prior$P0

    if(is.null(beta0)) beta0 <- rep(0, k)
    if(length(beta0) == 1) beta0 <- rep(beta0, k)
    if(is.null(B0)) B0 <- 100*diag(k)
    if(is.null(nuF)) nuF <- 4
    if(is.null(rho0)) rho0 <- q+4
    if(is.null(R0)) R0 <- 100*diag(q)
    if(is.null(nuG)) nuG <- 4
    if(is.null(nu00)) nu00 <- 4
    if(is.null(delta00)) delta00 <- 2
    if(is.null(nu0)) nu0 <- 4
    if(is.null(delta0)) delta0 <- 2
    if(is.null(phi0)) phi0 <- rep(0, ar+ma)
    if(is.null(P0)) P0 <- diag(ar+ma)
    if(ar==0&ma==0){
      phi0 <- 0
      P0 <- 1
      }

    ##check priors
    CheckNormalPriors(beta0, B0, k)
    CheckIGPriors(nu0, delta0)
    CheckWishartPriors(rho0, R0, q)
    if(hetero==TRUE) CheckGPriors(nu00, delta00)
    if(linkage=="t"){
      CheckGPriors(nuG, nuG)
      CheckGPriors(nuF, nuF)
      }
    if(ar!=0 |ma!=0) CheckNormalPriors(phi0, P0, ar+ma)

    ## check MCMC control parameters
    control <- CheckControl(control)
    burnin <- control$burnin
    iter <- control$iter
    iter2 <- control$iter2
    seed <- control$seed
    thin = control$thin
    verbose = control$verbose

    ##data supplements for C interface
    y <- as.double(y)
    X <- as.double(t(X))
    W <- as.double(t(W))
    dim <- as.integer(c(n, nn, k, q))
    ind <- as.integer(ind)
    beta0 <- as.double(beta0)
    B0 <- as.double(t(B0))
    nuF0 <- as.double(nuF)
    rho0 <- as.double(rho0)
    R0 <- as.double(t(R0))
    nuG0 <- as.double(nuG)
    nu00 <- as.double(nu00)
    delta00 <- as.double(delta00)
    nu0 <- as.double(nu0)
    delta0 <- as.double(delta0)
    phi0 <- as.double(phi0)
    P0 <- as.double(t(P0))
    mo <- floor(iter/thin)
    control <- as.integer(c(burnin, iter, mo, iter2, seed, thin, verbose))
    link <- as.integer(linkage=="t")
    hetero <- as.integer(hetero)
    arma <- as.integer(c(ar, ma))
    betam <- double(k*mo)
    bm <- double(mo*n*q)
    Dm <- double(mo*q*q)
    sigma2m <- double(mo)
    if(hetero==TRUE) sigma2m <- double(n*mo)
    phim <- double(1)
    if(ar!=0|ma!=0) phim <- double((ar+ma)*mo)
    lnmarglik <- double(1)

    result <- .C("bayespanel", control, link, hetero, arma,
                 y, X, W, ind,dim, beta0, B0, nuF0, rho0, R0, nuG0, nu00, delta00,
                 nu0, delta0, phi0, P0,
                 beta = betam, b = bm, D = Dm, sigma2 = sigma2m, phi = phim,
                 lnmarglik = lnmarglik)
    beta <- matrix(result$beta, nrow = mo, byrow = T)

    b <- matrix(result$b, nrow = mo, byrow = T)
    D <- matrix(result$D, nrow = mo, byrow = T)
    sigma2 <- matrix(result$sigma2, nrow = mo, byrow=T)
    if(ar!=0|ma!=0) phi <- matrix(result$phi, nrow = mo, byrow=T )
    beta <- mcmc(beta, start = burnin+1, end = burnin+iter, thin = thin)
    varnames(beta) <- as.list(Xnames)
    #b <- mcmc(b, start = burnin+1, end = burnin+iter, thin = thin)
    bmean <- apply(b, 2, mean)
    bmean <- matrix(bmean, nrow = n, byrow=T)
    dimnames(bmean)[[1]] <- as.list(indnames)
    dimnames(bmean)[[2]] <- as.list(Wnames)

    Dform <- Dform(Wnames)
    D <- mcmc(D[,Dform$Dind], start = burnin+1, end = burnin+iter,
    thin = thin)
    varnames(D) <- as.list(Dform$Dlabel)

    sigma2 <- mcmc(sigma2, start = burnin+1, end = burnin+iter, thin = thin)
    if(ar!=0|ma!=0) phi <- mcmc(phi, start = burnin+1, end = burnin+iter, thin = thin)

    cat("\n \n Summary Statistics for Posterior Distribution: \n")
    cat("\n \n The Model Specification is: \n")
    print(cl)
    cat("\n \n Posterior Estimates of Fixed Effect Coefficients:\n")
    print(summary(beta))
    cat("\n \n Posterior Mean Estimates of Random Coefficients:\n")
    print(bmean)
    cat("\n \n Posterior Estimates of Variance Structure for Random Coefficients:\n")
    print(summary(D))
    if (hetero==FALSE){
      cat("\n \n Posterior Estimates of Variance:\n")
      print(summary(sigma2))
    }

    out = list(Call = cl, beta = beta, b = bmean, D = D, sigma2 = sigma2)
    return(out)
}










