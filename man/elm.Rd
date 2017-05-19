\name{elm}
\alias{elm}

\title{ELM neural network.}

\description{Fit ELM neural network.}

\usage{
elm(y,hd=NULL,type=c("lasso","ridge","step","lm"),reps=20,
    comb=c("median","mean","mode"),lags=NULL,keep=NULL,difforder=NULL,
    outplot=c(FALSE,TRUE),sel.lag=c(TRUE,FALSE),
    direct=c(FALSE,TRUE),allow.det.season=c(TRUE,FALSE),
    det.type=c("auto","bin","trg"),xreg=NULL,xreg.lags=NULL,
    xreg.keep=NULL,barebone=c(FALSE,TRUE))
}
\arguments{
  \item{y}{
    Input time series. Can be ts or msts object. 
    }
  \item{hd}{
    Number of hidden nodes. This can be a vector. Use NULL to automatically specify.
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
  \item{lags}{
    Lags of y to use as inputs. If none provided then 1:frequency(y) is used. Use 0 for no univariate lags.
    }
    \item{keep}{
    Logical vector to force lags to stay in the model if sel.lag == TRUE. If NULL then it keep = rep(FALSE,length(lags)).
    }
    \item{difforder}{
    Vector including the differencing lags. For example c(1,12) will apply first and seasonal (12) differences. For no differencing use 0. For automatic selection use NULL.
    }
  \item{outplot}{
    Provide plot of model fit. Can be TRUE or FALSE.
    }
  \item{sel.lag}{
    Use selective lags only. Can be TRUE or FALSE. 
    }
  \item{direct}{
    Use direct input-output connections to model strictly linear effects. Can be TRUE or FALSE. 
    }
  \item{allow.det.season}{
    Permit modelling seasonality with deterministic dummies.
  }
  \item{det.type}{
    Type of deterministic seasonality dummies to use. This can be "bin" for binary or "trg" for a sine-cosine pair. With "auto" if ony a single seasonality is used and periodicity is up to 12 then "bin" is used, otherwise "trg".
  }
  \item{xreg}{
    Exogenous regressors. Each column is a different regressor and the sample size must be at least as long as the target in-sample set, but can be longer. 
  }
  \item{xreg.lags}{
    This is a list containing the lags for each exogenous variable. Each list is a numeric vector containing lags. If xreg has 3 columns then the xreg.lags list must contain three elements. If NULL then it is automatically specified. 
  }
  \item{xreg.keep}{
    List of logical vectors to force lags of xreg to stay in the model if sel.lag == TRUE. If NULL then all exogenous lags can be removed.
    }
  \item{barebone}{
    Use an alternative elm implementation (written in R) that is faster when the number of inputs is very high. Typically not needed.
  }
}
\value{An object of class "\code{elm}". If barebone == TRUE then the object inherits a second class "\code{elm.fast}"
The function \code{plot} produces a plot the network architecture.
An object of class \code{"elm"} is a list containing the following elements:
\item{net}{ELM networks. If it is of class "\code{elm.fast}" then this is NULL.}
\item{hd}{Number of hidden nodes. If it is of class "\code{elm.fast}" this is a vector with a different number for each repetition.}
\item{W.in}{NULL unless it is of class "\code{elm.fast}". Contains the input weights.}
\item{W}{Output layer weights for each repetition.}
\item{b}{Contains the output node bias for each training repetition.}
\item{W.dct}{Contains the direct connection weights if argument direct == TRUE. Otherwise is NULL.}
\item{lags}{Input lags used.}
\item{xreg.lags}{xreg lags used.}
\item{difforder}{Differencing used.}
\item{sdummy}{Use of deterministic seasonality.}
\item{ff}{Seasonal frequencies detected in data (taken from ts or msts object).}
\item{ff.det}{Seasonal frequencies coded using deterministic dummies.}
\item{det.type}{Type of determistic seasonality.}
\item{y}{Input time series.}
\item{minmax}{Scaling structure.}
\item{xreg.minmax}{Scaling structure for xreg variables.}
\item{comb}{Combination operator used.}
\item{type}{Estimation used for output layer weights.}
\item{direct}{Presence of direct input-output connections.}
\item{fitted}{Fitted values.}
\item{MSE}{In-sample Mean Squared Error.}
}
\references{
    \itemize{
        \item{For combination operators see: Kourentzes N., Barrow B.K., Crone S.F. (2014) \href{http://kourentzes.com/forecasting/2014/04/19/neural-network-ensemble-operators-for-time-series-forecasting/}{Neural network ensemble operators for time series forecasting}. \emph{Expert Systems with Applications}, \bold{41}(\bold{9}), 4235-4244.}
        \item{For ELMs see: Huang G.B., Zhou H., Ding X. (2006) Extreme learning machine: theory and applications. \emph{Neurocomputing}, \bold{70}(\bold{1}), 489-501.}
}
}
\note{To use elm with Temporal Hierarchies (\href{https://cran.r-project.org/package=thief}{thief} package) see \code{\link{elm.thief}}.
The elm function by default calls the \code{neuralnet} function. If barebone == TRUE then it uses an alternative implementation (\code{TStools:::elm.fast}) which is more appropriate when the number of inputs is several hundreds.}
\author{
Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}
}
\seealso{
\code{\link{forecast.elm}}, \code{\link{plot.elm}}, \code{\link{elm.thief}}, \code{\link{mlp}}.
}
\examples{
fit <- elm(AirPassengers)
print(fit)
plot(fit)
frc <- forecast(fit,h=36)
plot(frc)
}
\keyword{ ~elm }
\keyword{ ~thief }
\keyword{ ~ts }