\name{mlp}
\alias{mlp}

\title{Fit MLP neural network.}

\description{Fit MLP neural network.}

\usage{
mlp(y,m=frequency(y),hd=NULL,reps=20,
    comb=c("median","mean","mode"),lags=NULL,
    difforder=-1,outplot=c(FALSE,TRUE),
    sel.lag=c(TRUE,FALSE),allow.det.season=c(TRUE,FALSE),
    det.type=c("auto","bin","trg"),
    xreg=NULL,xreg.lags=NULL,
    hd.auto.type=c("set","valid","cv","elm"),...)
}

\arguments{
  \item{y}{
    Input time series. Can be ts or msts object.
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
  \item{hd.auto.type}{
    Used only if hd==NULL. "set" fixes hd=5. "valid" uses a 20\% validation set (randomly) sampled to find the best number of hidden nodes. "cv" uses 5-fold cross-validation. "elm" uses ELM to estimate the number of hidden nodes (experimental).
  }
  \item{...}{
    Additional inputs for neuralnet function. 
    }
}
\value{
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
}
\references{
    \itemize{
        \item{For combination operators see: Kourentzes N., Barrow B.K., Crone S.F. (2014) Neural network ensemble operators for time series forecasting. \emph{Expert Systems with Applications}, \bold{41}(\bold{9}), 4235-4244.}
}
}
\author{
Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}
}
\examples{
y <- log(AirPassengers)
fit <- mlp(y)
frc <- forecast(fit,h=36,outplot=TRUE)
}
