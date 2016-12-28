\name{elm}
\alias{elm}

\title{Fit ELM neural network.}

\description{Fit ELM neural network.}

\usage{
elm(y,hd=50,type=c("lasso","step","lm"),reps=20,
    comb=c("median","mean","mode"),lags=NULL,
    difforder=-1,outplot=c(FALSE,TRUE),sel.lag=c(FALSE,TRUE),
    direct=c(FALSE,TRUE))
}

\arguments{
  \item{y}{
    Input time series. 
    }
  \item{hd}{
    Number of hidden nodes. This can be a vector.
    }
  \item{type}{
    Estimation type for output layer weights. Can be "lasso" (lasso with CV), "step" (stepwise regression with AIC) or "lm" (linear regression).
    }
   \item{reps}{
    Number of networks to train.
    }
  \item{comb}{
    Combination operator for forecasts when reps > 1. Can be "median", "mode" (based on KDE estimation) and "mean".
    }
  \item{lags}{
    Lags of y to use as inputs. If none provided then 1:frequency(y) is used.
    }
    \item{difforder}{
    Vector including the differencing lags. For example c(1,12) will apply first and seasonal (12) differences. For no differencing use NULL. For automatic selection use -1.
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
}
\value{
\item{net}{ELM networks.}
\item{hd}{Number of hidden nodes.}
\item{W}{Output layer weights for ech repetition.}
\item{lags}{Input lags used.}
\item{difforder}{Differencing used.}
\item{y}{Input time series.}
\item{minmax}{Scaling structure.}
\item{comb}{Combination operator used.}
\item{type}{Estimation used for output layer weights.}
\item{direct}{Presence of direct input-output connections.}
\item{fitted}{Fitted values.}
\item{MSE}{In-sample Mean Squared Error.}
}
\references{
    \itemize{
        \item{For combination operators see: Kourentzes N., Barrow B.K., Crone S.F. (2014) Neural network ensemble operators for time series forecasting. \emph{Expert Systems with Applications}, \bold{41}(\bold{9}), 4235-4244.}
        \item{For ELMs see: ... .}
}
}
\author{
Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}
}
\examples{
y <- log(AirPassengers)
fit <- elm(y)
frc <- forecast(fit,h=36,outplot=TRUE)
}
