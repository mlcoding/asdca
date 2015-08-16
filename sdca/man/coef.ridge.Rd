\name{coef.ridge}
\alias{coef.ridge}
\title{
Extract Model Coefficients for an object with S3 class \code{"ridge"}
}
\description{
Extract estimated regression coefficient vectors from the solution path.
}
\usage{
\method{coef}{ridge}(object, lambda.idx = c(1:3), beta.idx = c(1:3), \dots)
}
\arguments{
  \item{object}{
  An object with S3 class \code{"ridge"} 
  }
  \item{lambda.idx}{
  The indices of the regularizaiton parameters in the solution path to be displayed. The default values are \code{c(1:3)}.
  }
  \item{beta.idx}{
  The indices of the estimate regression coefficient vectors in the solution path to be displayed. The default values are \code{c(1:3)}.
  }
  \item{\dots}{
  Arguments to be passed to methods.
  }
}

\author{
Xingguo Li, Tuo Zhao, Daniel Vainsencher and Han Liu \cr
Maintainer: Xingguo Li <xingguo.leo@gmail.com>
}

\seealso{
\code{\link{stocml}} and \code{\link{stocml-package}}.
}
