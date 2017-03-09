\name{elm.fast}
\alias{elm.fast}

\title{ELM (fast) neural network.}

\description{Fit ELM (fast) neural network. This is an ELM implementation that does not rely on neuralnets package.}

\usage{
elm.fast(y,x,hd=NULL,type=c("lasso","ridge","step","ls"),reps=20,
                     comb=c("median","mean","mode"),direct=c(FALSE,TRUE),
                     linscale=c(TRUE,FALSE),output=c("linear","logistic"),
                     core=c("FALSE","TRUE"),ortho=c(FALSE,TRUE))
}
\arguments{
  \item{y}{
    Target variable.
    }
  \item{x}{
    Explanatory variables. Each column is a variable. 
  }
  \item{hd}{
    Starting number of hidden nodes (scalar). Use NULL to automatically specify.
    }
  \item{type}{
    Estimation type for output layer weights. Can be "lasso" (lasso with CV), "ridge" (ridge regression with CV), "step" (stepwise regression with AIC) or "lm" (linear regression).
    }
   \item{reps}{
    Number of networks to train.
    }
  \item{comb}{
    Combination operator for forecasts when reps > 1. Can be "median", "mode" (based on KDE estimation) and "mean".
    }
  \item{direct}{
    Use direct input-output connections to model strictly linear effects. Can be TRUE or FALSE. 
    }
  \item{linscale}{
    Scale inputs linearly between -0.8 to 0.8. If output == "logistic" then scaling is between 0 and 1.
  }
  \item{output}{
    Type of output layer. It can be "linear" or "logistic". If "logistic" then type must be set to "lasso".
  }
  \item{core}{
    If TRUE skips calculation of final fitted values and MSE. Called internally by \code{"elm"} function.
  }
  \item{ortho}{
    If TRUE then the initial weights between the input and hidden layers are orthogonal (only when number of input variable <= sample size).
  }
}
\value{An object of class "\code{elm.fast}".
The function \code{plot} produces a plot the network fit.
An object of class \code{"elm.fast"} is a list containing the following elements:
\item{hd}{Number of hidden nodes. This is a vector with a different number for each training repetition.}
\item{W.in}{Input weights for each training repetition.}
\item{W}{Output layer weights for each training repetition.}
\item{b}{Output node bias for each training repetition.}
\item{W.dct}{Direct connection weights argument if direct == TRUE for each training repetition. Otherwuse NULL.}
\item{fitted.all}{Fitted values for each training repetition.}
\item{fitted}{Ensemble fitted values.}
\item{y}{Target variable.}
\item{type}{Estimation used for output layer weights.}
\item{comb}{Combination operator used.}
\item{direct}{Presence of direct input-output connections.}
\item{output}{Type of output node.}
\item{minmax}{If scaling is used this contains the scaling information for the target variable.}
\item{minmax.x}{If scaling is used this contains the scaling information for the input variables.}
\item{MSE}{In-sample Mean Squared Error.}
}
\references{
    \itemize{
        \item{For combination operators see: Kourentzes N., Barrow B.K., Crone S.F. (2014) \href{http://kourentzes.com/forecasting/2014/04/19/neural-network-ensemble-operators-for-time-series-forecasting/}{Neural network ensemble operators for time series forecasting}. \emph{Expert Systems with Applications}, \bold{41}(\bold{9}), 4235-4244.}
        \item{For ELMs see: Huang G.B., Zhou H., Ding X. (2006) Extreme learning machine: theory and applications. \emph{Neurocomputing}, \bold{70}(\bold{1}), 489-501.}
}
}
\note{This implementation of ELM is more appropriate when the number of inputs is several hundreds. For time series modelling use \code{\link{elm}} instead.}
\author{
Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}
}
\seealso{
\code{\link{predict.elm.fast}}, \code{\link{plot.elm.fast}}, \code{\link{elm}}.
}
\examples{
p <- 2000
n <- 150
X <- matrix(rnorm(p*n),nrow=n)
b <- cbind(rnorm(p))
Y <- X \%*\% b
fit <- elm.fast(Y,X)
print(fit)
}
\keyword{ ~elm }