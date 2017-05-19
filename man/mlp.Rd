\name{mlp}
\alias{mlp}

\title{MLP neural network.}

\description{Fit MLP neural network.}

\usage{
mlp(y,m=frequency(y),hd=NULL,reps=20,comb=c("median","mean","mode"),
    lags=NULL,keep=NULL,difforder=NULL,outplot=c(FALSE,TRUE),
    sel.lag=c(TRUE,FALSE),allow.det.season=c(TRUE,FALSE),
    det.type=c("auto","bin","trg"),xreg=NULL, xreg.lags=NULL,
    xreg.keep=NULL,hd.auto.type=c("set","valid","cv","elm"),
    hd.max=NULL, ...)
}

\arguments{
  \item{y}{
    Input time series. Can be ts or msts object.
    }
  \item{m}{
    Frequency of the data.
    }
  \item{hd}{
    Number of hidden nodes. This can be a vector.
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
  \item{hd.auto.type}{
    Used only if hd==NULL. "set" fixes hd=5. "valid" uses a 20\% validation set (randomly) sampled to find the best number of hidden nodes. "cv" uses 5-fold cross-validation. "elm" uses ELM to estimate the number of hidden nodes (experimental).
  }
  \item{hd.max}{
    When hd.auto.type is set to either "valid" or "cv" then this argument can be used to set the maximum number of hidden nodes to evaluate, otherwise the maximum is set automatically.
    }
  \item{...}{
    Additional inputs for neuralnet function. 
    }
}
\value{An object of class "\code{mlp}".
The function \code{plot} produces a plot the network architecture.
An object of class \code{"mlp"} is a list containing the following elements:
\item{net}{MLP networks.}
\item{hd}{Number of hidden nodes.}
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
\item{fitted}{Fitted values.}
\item{MSE}{In-sample Mean Squared Error.}
\item{MSEH}{If hd.auto.type is set to either "valid" or "cv" an array of the MSE error for each network size is provided. Otherwise this is NULL.}
}
\references{
    \itemize{
        \item{For combination operators see: Kourentzes N., Barrow B.K., Crone S.F. (2014) \href{http://kourentzes.com/forecasting/2014/04/19/neural-network-ensemble-operators-for-time-series-forecasting/}{Neural network ensemble operators for time series forecasting}. \emph{Expert Systems with Applications}, \bold{41}(\bold{9}), 4235-4244.}
}
}
\note{To use mlp with Temporal Hierarchies (\href{https://cran.r-project.org/package=thief}{thief} package) see \code{\link{mlp.thief}}.}
\author{
Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}
}
\seealso{
\code{\link{forecast.mlp}}, \code{\link{plot.mlp}}, \code{\link{mlp.thief}}, \code{\link{elm}}.
}
\examples{
fit <- mlp(AirPassengers)
print(fit)
plot(fit)
frc <- forecast(fit,h=36)
plot(frc)
}
\keyword{ ~mlp }
\keyword{ ~thief }
\keyword{ ~ts }
