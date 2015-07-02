\name{stocml}
\alias{stocml}

\title{
Stochastic Machine Learning
}

\description{
The function "stocml" implements The user interface.
}

\usage{
stocml(X, Y, lambda = NULL, nlambda = NULL, lambda.min.ratio = NULL,
        family = "gaussian", alg = "sdca", standardize = FALSE, design.sd = TRUE,
        res.sd = FALSE, prec = 1e-4, max.ite = 1e4, verbose = TRUE)
}

\arguments{
  \item{X}{For sparse linear regression and sparse logistic regression, \code{X} is an \eqn{n} by \eqn{d} design matrix. For sparse multivariate regression, \code{X} is an \eqn{n} by \eqn{d} design matrix. For sparse inverse column operation, there are 2 options when \code{family = "npn"}: (1) \code{data} is an \code{n} by \code{d} data matrix (2) a \code{d} by \code{d} sample covariance matrix. The program automatically identifies the input matrix by checking the symmetry. (\code{n} is the sample size and \code{d} is the dimension). }
  \item{Y}{For sparse linear regression and sparse logistic regression, \code{Y} is an \eqn{n} dimensional response vector. For sparse multivariate regression, \code{Y} is an \eqn{n} by \eqn{m} response matrix. For sparse inverse column operation, no input for \code{Y} is needed.}
  \item{lambda}{A sequence of decresing positive value to control the regularization. Typical usage is to leave the input \code{lambda = NULL} and have the program compute its own \code{lambda} sequence based on \code{nlambda} and \code{lambda.min.ratio}. Users can also specify a sequence to override this. Default value is from \eqn{lambda.max} to \code{lambda.min.ratio*lambda.max}. For Lq regression, the default value of \eqn{lambda.max} is \eqn{\pi\sqrt{\log(d)/n}}. The default value of \eqn{lambda.max} is the minimum regularization parameter which yields an all-zero estiamtes.}
  \item{nlambda}{The number of values used in \code{lambda}. Default value is 5.}
  \item{lambda.min.ratio}{The smallest value for \code{lambda}, as a fraction of the uppperbound (\code{MAX}) of the regularization parameter. The program can automatically generate \code{lambda} as a sequence of length = \code{nlambda} starting from \code{MAX} to \code{lambda.min.ratio*MAX} in log scale. The default value is \code{0.25}.}
  \item{family}{Options for model. Sparse linear regression and sparse multivariate regression is applied if \code{family = "gaussian"}, sparse logistic regression is applied if \code{family = "logit"} and sparse column inverse operation is applied if \code{family = "npn"}. The default value is \code{"gaussian"}.}
  \item{alg}{Options for active set identification. Cyclic search is applied if \code{alg = "cyclic"}, greedy search is applied if \code{alg = "greedy"}, proximal gradient based search is applied if \code{alg = "prox"} and stochastic search is applied if \code{alg = "stoc"}. The default value is \code{"cyclic"}.}
  \item{standardize}{Variables are standardized to have mean zero and unit standard deviation if \code{standardize = TRUE}. The default value is \code{FALSE}.}
  \item{design.sd}{Flag of whether the design matrix is standardized. The default value is \code{TRUE}.}
  \item{res.sd}{Flag of whether the response varialbles are standardized. The default value is \code{FALSE}.}
  \item{prec}{Stopping criterion. The default value is 1e-4.}
  \item{max.ite}{The iteration limit. The default value is 1e4.}
  \item{verbose}{Tracing information is disabled if \code{verbose = FALSE}. The default value is \code{TRUE}.}
  }
\details{

}

\value{

}

\author{
Xingguo Li, Tuo Zhao and Daniel Vainsencher \cr
Maintainer: Xingguo Li <xingguo.leo@gmail.com>
}

\references{

}

\seealso{
\code{\link{stocml-package}}.
}

\examples{
}