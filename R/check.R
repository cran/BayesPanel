################################################################################
## Functions for checking arguments in the model.
## Not accessible by the end user
## Generic for every function
################################################################################

"CheckIndex" <- function(data, index){
    # index can take two character argument of names for subject and time
    if(is.null(index)) stop("Please specify the index argument and recall.\n")
    if(match(index[1], colnames(data),0) == 0){
         stop("The subject index variable", index[1], "can not be found.\n")
         }
    if(length(index) == 2 & match(index[2], colnames(data), 0) == 0){
         stop ("The time index varialbe", index[2], "can not be found.\n")
         }
    sind <- match(index[1], colnames(data))
    order <- order(data[,sind])
    if (length(index) == 2){
      tind <- match(index[2], colnames(data))
      order <- order(data[,sind], data[,tind])
    }
    if(max(order - 1:nrow(data)) > 0 | min(order - 1:nrow(data)) <0){
      warning("The data frame is not ordered according to the index,
               sorting was done automatically.\n")
      }
    data <- data[order,]
    ind <- as.numeric(table(data[,sind]))
    indnames <- paste(unique(data[,sind]))
    return(list(data = data, ind = ind, indnames = indnames ))
           ## can also return order, and sorting back after estimation
  }

"CheckARMA" <- function(arma){
  if(length(arma)>2) warning ("More than two values for arma argument is supplemented, only the first two values will be used.\n")
  ar <- arma[1]
  ma <- arma[2]

  if (ar>0 | ma>0) warning("Models with time series error structure are notimplemented in the current version, the model will be estimated with independent error structure")
 #  if (ar>2 | ma>2) stop
#  ("Models with error structure beyond ar(2)ma(2) are not implemented, please respecify.\n")
  if(ar<0 | ma<0 | ar-floor(ar)>0 | ma-floor(ma)>0) stop("Argumetns for arma must be positive integers, please respecify.\n")
  return(c(ar, ma))
  }


"CheckControl" <- function(control){
   if (is.null(control$seed)) control$seed = 1234567
   if (is.null(control$burnin)) control$burnin = 1000
   if (is.null(control$iter)) control$iter = 10000
   if (is.null(control$iter2)) control$iter2 = 10000
   if (is.null(control$thin)) control$thin =1
   if (is.null(control$verbose)) control$verbose = 1000
   if (control$iter %% control$thin != 0){
     warning ("Number of MCMC iterations is not divisible by thinning
          interval.\n")
     }
   if (control$seed %% 1 != 0 | control$seed < 0){
     stop(" The seed argument must be a non-negative integer,\n
       Please respecify the value of \"burnin\" and run again.\n")
     }
   if (control$burnin %% 1 != 0 | control$burnin < 0){
     stop(" Number of burn-in iterations must be a non-negative integer,\n
       Please respecify the value of \"burnin\" and run again.\n")
     }
   if (control$iter %% 1 != 0 | control$iter < 0 ){
     stop("Number of MCMC iterations must be a non-negative integer,\n
       Please respecify the value of \"iter\" and run again.\n")
     }
   if (control$iter2 %% 1 != 0 | control$iter2 < 0 ){
     stop("Number of additional MCMC iterations for marginal likelihood calculation,\n
       must be a non-negative integer,\n
       Please respecify the value of \"iter2\" and run again.\n")
     }
   if (control$thin %% 1 != 0 | control$thin <= 0 ){
     stop("Value of thinning interval must be a positive integer,\n
       Please respecify the value of \"thin\" and run again.\n")
     }
   if (control$verbose %% 1 != 0 | control$verbose < 0 ){
     stop("Value of verbose interval must be a positive integer,\n
       Please respecify the value of \"verbose\" and run again.\n")
     }
   return(control)
   }

"CheckNormalPriors" <- function(mu, Sigma, dim){
    if (length(mu) != dim){
      stop("The prior mean ", deparse(substitute(mu)),
        " is not in the right dimension of ", dim , ", \n",
        "Please respecify ", deparse(substitute(mu)), " and run again.\n")
      }
    if (dim(Sigma)[1] != dim(Sigma)[2]){
      stop("The prior variance ", deparse(substitute(Sigma)),
        " must be a squared matrix,\n",
        "Please respecify ", deparse(substitute(Sigma)), " and run again.\n")
      }
    if (dim(Sigma)[1] != dim){
      stop("The prior variance ", deparse(substitute(Sigma)),
        " is not in the right dimension of ", dim, ",\n",
        "Please respecify ", deparse(substitute(Sigma)), " and run again.\n")
      }
    if (min(eigen(Sigma)$values) < 0){
      stop("The prior variance ", deparse(substitute(Sigma)),
        " is not positive definite,\n",
        "Please respecify ", deparse(substitute(Sigma)), " and run again.\n")
      }
    return(0)
    }

"CheckIGPriors" <- function(nu, delta){
    if (nu < 0){
      stop("The prior parameter ", deparse(substitute(nu)),
      " must be positive,\n", "Please respecify ",
      deparse(substitute(nu)), " and run again.\n")
      }
    if (delta < 0){
      stop("The prior parameter ", deparse(substitute(delta)),
      " must be positive,\n", "Please respecify ",
      deparse(substitute(delta)), " and run again.\n")
      }
    return (0)
    }


"CheckGPriors" <- function(nu, delta){
    if (nu < 0){
      stop("The prior parameter ", deparse(substitute(nu)),
      " must be positive,\n", "Please respecify ",
      deparse(substitute(nu)), " and run again.\n")
      }
    if (delta < 0){
      stop("The prior parameter ", deparse(substitute(delta)),
      " must be positive,\n", "Please respecify ",
      deparse(substitute(delta)), " and run again.\n")
      }
    return (0)
    }

"CheckWishartPriors" <- function(rho, R, dim){
    if (rho < dim+1){
      stop("The prior parameter ", deparse(substitute(rho)),
      " must be a real value great than", deparse(substitute(dim)), "+1\n", "Please respecify ", deparse(substitute(rho)), " and run again.\n")
    }
    if (dim(R)[1] != dim(R)[2]){
      stop("The prior variance ", deparse(substitute(R)),
        " must be a squared matrix,\n",
        "Please respecify ", deparse(substitute(R)), " and run again.\n")
      }
    if (dim(R)[1] != dim){
      stop("The prior variance ", deparse(substitute(R)),
        " is not in the right dimension of ", dim, ",\n",
        "Please respecify ", deparse(substitute(R)), " and run again.\n")
      }
    if (min(eigen(R)$values) < 0){
      stop("The prior variance ", deparse(substitute(R)),
        " is not positive definite,\n",
        "Please respecify ", deparse(substitute(R)), " and run again.\n")
      }
    return(0)
    }

"Dform" <- function(Wnames){
  q <- length(Wnames)
  Dind <- matrix(1:q^2, q, q)
  Dind <- Dind[lower.tri(Dind, diag=TRUE)]
  Dlabel1 <- rep(1:q, 1:q)
  Dlabel1 <- Wnames[Dlabel1]
  Dlabel2 <- NULL
  for (i in 1:q){
    Dlabel2 <- c(Dlabel2, 1:i)
  }
  Dlabel2 <- Wnames[Dlabel2]
  Dlabel <- paste(Dlabel1, Dlabel2, sep=":")
  return(list(Dind = Dind, Dlabel= Dlabel))
}
