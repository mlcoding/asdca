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
       method = "ridge", alg = "sdca", design.sd = TRUE,
       res.sd = FALSE, prec = 1e-4, max.ite = 1e4, verbose = TRUE)
}

\arguments{
  \item{X}{An \eqn{n} by \eqn{d} design matrix. (\code{n} is the sample size and \code{d} is the dimension). }
  \item{Y}{For linear regression, \code{Y} is an \eqn{n} dimensional response vector. For linear classification, \code{Y} is an \eqn{n} dimensional 1 or -1 category vector.}
  \item{lambda}{A sequence of decresing positive values to control the penalty. Typical usage is to leave the input \code{lambda = NULL} and have the program compute its own \code{lambda} sequence based on \code{nlambda} and \code{lambda.min.ratio}. Users can also specify a sequence to override this. Default value is from \eqn{lambda.max} to \code{lambda.min.ratio*lambda.max}.}
  \item{nlambda}{The number of values used in \code{lambda}. Default value is 5.}
  \item{lambda.min.ratio}{The smallest value for \code{lambda}, as a fraction of the uppperbound (\code{MAX}) of the regularization parameter. The program can automatically generate \code{lambda} as a sequence of length = \code{nlambda} starting from \code{MAX} to \code{lambda.min.ratio*MAX} in log scale. The default value is \code{0.25}.}
  \item{method}{Options for the method. Linear regression is applied if \code{method = "ridge"} and linear classification is applied if \code{method = "svm"}. The default value is \code{"ridge"}.}
  \item{alg}{Options for active set identification. The default value is \code{"sdca"}.}
  \item{design.sd}{Flag of whether the design matrix is standardized for linear regression. The default value is \code{TRUE}.}
  \item{res.sd}{Flag of whether the response varialbles are standardized for linear regression. The default value is \code{FALSE}.}
  \item{prec}{Stopping criterion. The default value is 1e-4.}
  \item{max.ite}{The iteration limit. The default value is 1e4.}
  \item{verbose}{Tracing information is disabled if \code{verbose = FALSE}. The default value is \code{TRUE}.}
  }
\details{

}

\value{
An object with S3 classes \code{"ridge"} and \code{"class"} corresponding to linear regression and linear classification is returned: 
  \item{beta}{
  A matrix of regression estimates whose columns correspond to regularization parameters.
  }
  \item{intercept}{The value of intercepts corresponding to regularization parameters for linear regression. }
  \item{Y}{The value of \code{Y} used in the program. }
  \item{X}{The value of \code{X} used in the program. }
  \item{lambda}{
  The sequence of regularization parameters \code{lambda} used in the program.
  }
  \item{nlambda}{
  The number of values used in \code{lambda}.
  }
  \item{method}{
  The \code{method} from the input. 
  }
  \item{method}{
  The \code{method} from the input. 
  }
  \item{alg}{
  The \code{alg} from the input. 
  }
  \item{ite}{
  A list of vectors where the i-th entries of ite[[1]] and ite[[2]] correspond to the outer iteration and inner iteration of i-th regularization parameter respectively.
  }
  \item{verbose}{The \code{verbose} from the input.}
}

\author{
Xingguo Li, Tuo Zhao, Daniel Vainsencher and Han Liu \cr
Maintainer: Xingguo Li <xingguo.leo@gmail.com>
}

\references{
1. Shai Shalev-Shwartz and Tong Zhang. Accelerated Proximal Stochastic Dual Coordinate Ascent for Regularized Loss Minimization.
}

\seealso{
\code{\link{stocml-package}}.
}

\examples{
n = 500
d = 20
X = matrix(rnorm(n*d), n, d)
beta = rnorm(d)
lambda = c(1e-2,1e-3,1e-4)
Y = X\%*\%beta + rnorm(n)
out = stocml(X, Y, lambda=lambda)
}