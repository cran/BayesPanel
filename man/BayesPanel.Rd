\name{BayesPanel}
\alias{BayesPanel}
\title{Bayesian Methods for Random Effects Panel Data Regression}
\description{
  This function generates posterior distribution of parameters for a
  random effects panel data regression model. 
  }
  
\usage{
BayesPanel(formula, data, subset, na.action, index,
           linkage = c("normal", "t"), hetero = FALSE, arma = c(0, 0),
           prior = list(beta0 = NULL, B0 = NULL, nuF = NULL,
                        rho0 = NULL, R0 = NULL,
                        nuG = NULL, nu00 = NULL, delta00 = NULL,
                        nu0 = NULL,delta0 = NULL, phi0 = NULL, P0 = NULL),
           control = list(seed = 314159265, burnin = 1000, iter = 10000,
                          iter2 = 10000, thin = 1, verbose = 0),
           ...  )
}

\arguments{
    \item{formula}{Model formula. We use an extended expression of
    formula in the Formula package. The response varialbe is in the left
    hand side of "~", followed by the fixed effects varaibles X and random
    effects varaibles W in the right hand side. X and W are further
    separated by "|". For example, a valid formula would look like  y ~ X1 + X2 | W1 + W2 .}

  \item{data}{The data frame used.}
  
  \item{subset}{The subset to use for the data frame.}

  \item{na.action}{Method to be used for na values. Not implemented in
  the current version.}

  \item{index}{Strings of the variable names for subject and time. If
    supplemented with one string, then it is treated as the variable
    name for subject. For example, index = c("subject", "time").}

  \item{linkage}{The linkage distribution for the error term, either
    "normal" for normal distribution or "t" for student t distribution.}

  \item{hetero}{TRUE if the error term for each individual is
    heterogeneous, FALSE if they are homogeneous}

  \item{arma}{Implementation for time series dependence for the error
    term. Specify the options for arma model inference. First value is
    for ar and second value is for ma, eg. arma = c(1, 0) means ar(1)
    model. In this version, time series error structure has not been
    implemented.}

  \item{prior}{A list for the values of priors, see details.}

  \itemize{
    \item beta0, B0: priors for \eqn{\beta}.

    \item nuF: priors for \eqn{\eta}.
    
    \item rho0, R0: priors for \eqn{D}.

    \item nuG: priors for \eqn{\lambda}

    \item nu00, delta00: priors for \eqn{\delta_{00}}.

    \item nu0, delta0: priors for \eqn{\sigma^{2}}.

    \item phi0, P0: priors for arma coefficients.

    
  }

  \item{control}{A list of control variable for the MCMC sampling
    process.}

  \itemize{
    \item seed: seed used for random number generation.
    
    \item burnin: the number of burnin iterations to be used.

    \item iter: the number of iterations after burnin to be used.

    \item thin: the thinning interval to be used in the
    simulation. 
    
    \item verbose: a switch to determinies whether the progress should
  be printed onto screen. Default value of 0 means no printing. If
  \code{verbose} is greater than 0, then the process will be printed
  every \code{verbose}th iteration.
  }

  \item{...}{The further arguments to be passed.}
}
  
\value{
   An mcmc object that contains the posterior sample.  This 
   object can be summarized by functions provided by the coda package.
}

\details{
\code{BayesPanel} uses MCMC techniques to estimate parameters of a
random effect panel regression model.

The model specification is as follows:

\deqn{y_{i} = X_{i}\beta + W_{i}b_{i} + \epsilon_{i}}
\deqn{\epsilon_{i}|\sigma_{i}^{2}, \lambda_{i} \sim N_{n_{i}}(0, \sigma_{i}^{2}\lambda_{i}^{-1}\Omega_{i})}
\deqn{b_{i}|\eta_{i},D \sim N_{q} (0, \eta_{i}^{-1}D), i \le N}
\deqn{\beta\sim N_{k}(\beta_{0}, B_{0})}

  
}
\references{
  Chib, S. and B. P. Carlin. 1999. On MCMC sampling in hierarchical
  longitudinal models. \emph{Statistics and Computing.} 9, 17-26.

  Chib, S. 2008. Panel Data Modeling and Inference: A Bayesian
  Primer. \emph{The Econometrics of Panel Data. } (eds L. Matyas and
  P. Sevestre). 479-515. Springer-Verlag, Berlin Heidelberg. 
  
  Wakefield, J. C., A. F. M. Smith, A. Racine Poon, and A.E. Gelfand
  1994. Bayesian analysis of linear and non-linear population models by
  using the Gibbs sampler. \emph{Applied Statistics.} 43, 201-221.

  
}

\author{
  Chunhua Wu \email{chunhuawu@wustl.edu}
  
  Siddhartha Chib \email{chib@wustl.edu}
  }


\examples{
   \dontrun{
   library(plm)
   data(Gasoline)
   est <- BayesPanel(lgaspcar~lincomep+lrpmg+lcarpcap|lincomep,
                     data=Gasoline, index = c("country"),
                     control = list(verbose=10, thin=100) )
   }
}

\keyword{models}

\seealso{\code{\link[coda]{plot.mcmc}},\code{\link[coda]{summary.mcmc}}}

