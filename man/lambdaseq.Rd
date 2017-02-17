\name{lambdaseq}
\alias{lambdaseq}

\title{
Lambda Sequence generator for LASSO estimation
}

\description{
Calculates the lambdaMax value, which is the penalty term (lambda) beyond which coefficients are guaranteed to be all zero using 0*OLS solution and provides sequence of nLambda values to lambdaMin in logarithmic descent.
}


\usage{
lambdaseq(x, y, weight = NA, alpha = 1, standardise = T,
          lambdaRatio = 0.0001, nLambda = 100, addZeroLambda = F)
}

\arguments{
  \item{x}{
  Matrix of x variables
  }
  \item{y}{
  Vector of y values
  }
  \item{weight}{
  Vector of length(nrow(y)) for weighted LASSO estimation
  }
  \item{alpha}{
  Elastic net mixing value (default 1, range (0,1])
  }
  \item{standardise}{
  Standardises variables
  }
  \item{lambdaRatio}{
  Defines the ratio between lambdaMax and lambdaMin
  }
  \item{nLambda}{
  Defines lenght of lambda sequence
  }
  \item{addZeroLambda}{
  Sets last observation in lambda sequence to 0 which is the OLS solution
  }
}
\details{
Parameter \code{alpha} allows to define the elastic net mixing value or the relative balance between L2 and L1 penalty. Alpha = 1 --> lasso, otherwise elastic net. Alpha near zero --> nearly ridge regression.
}
\value{

  \item{lambda}{
  Sequence of lambda values from \code{lambdaMax} to \code{lambdaMin}
  }
  \item{lambdaMin}{
  Minimal lambda value
  }
  \item{lambdaMax}{
  Maximal lambda value
  }
  \item{nullMSE}{
  Mean Squared Error (MSE) of the fit
  }
}
\references{
Hastie, T., Tibshirani, R. and Wainwright, M. 2015. Statistical Learning with Sparsity. The Lasso and Generalizations. Boca Raton: CRC Press.
}
\author{
Oliver Schaer, \email{info@oliverschaer.ch},
Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}.
}
\note{
\code{nullMSE} is the MSE of the fit using just a constant term. It is provided in this function as a convenience, because it needs to be calculated in the same context as \code{lambdaMax} whenever \code{lambdaMax} is calculated.
}


\examples{
# Get some data from mtcars
y <- mtcars[,1]
x <- as.matrix(mtcars[,2:11])

# Create lambda sequence
lambda <- lambdaseq(x, y, standardise = T, lambdaRatio = 0.0001,
                    nLambda = 100, addZeroLambda = F)$lambda

# fit and obtain coefficients from LASSO model using glmnet package
library(glmnet)

fit.lasso <- cv.glmnet(y = y, x = x, family = "gaussian", alpha = 1,
                       standardize = T, lambda = lambda, intercept = T)

coef.lasso <- coef(fit.lasso, s = "lambda.1se")
}
\keyword{ ~lasso }

