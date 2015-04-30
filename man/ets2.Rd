\name{ets2}
\alias{ets2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Exponential smoothing in SSOE state-space model
}
\description{
Function constructs ETS model and returns forecast, fitted values, errors and matrix of states.
}
\usage{
ets2(data, model="ZZZ", persistence=NULL, phi=NULL,
    bounds=c("usual","admissible"),
    initial=NULL, initial.season=NULL,
    IC=c("AIC","AICc","BIC"), trace=FALSE, CF.type=c("GV","TLV","TV"),
    intervals=FALSE, int.w=0.95, xreg=NULL,
    holdout=FALSE, h=10, silent=FALSE, legend=TRUE,
    ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ the vector or ts object, containing the data needed to be forecasted.
  }
  \item{model}{ the type of ETS model. Can consist of 3 or 4 chars: ANN, AAN, AAdN, AAA, AAdA, MAdM etc.
  }
  \item{persistence}{ The vector of smoothing parameters. If NULL, the parameters are estimated.
  }
  \item{phi}{ The value of damping parameter. If NULL then the parameter is estimated.
  }
  \item{bounds}{ What type of bouds to use for the smoothing parameters. ONLY USUAL BOUNDS ARE IMPLEMENTED FOR NOW.
  }
  \item{initial}{ The vector of initial values for level and trend. If NULL, it is estimated in the optimization.
  }
  \item{initial.season}{ The vector of initial values for seasonal components. If NULL, they are estimated in the optimization.
  }
  \item{IC}{ The information criterion to use in the model selection. NOT IMPLEMENTED YET
  }
  \item{trace}{ If TRUE, the trace likelihood is used in the estimation with trace of the length h.
  }
  \item{CF.type}{ Type of Cost Function used in trace optimization. Can be "GV" - Generalized Variance, "TLV" - Total Logarythmic Variation and "TV" - Total Variation.
  }
  \item{intervals}{ If TRUE, the prediction interval is constructed. NOT IMPLEMENTED YET
  }
  \item{int.w}{ Defines the width of the prediction interval. NOT IMPLEMENTED YET
  }
  \item{xreg}{ The vector / matrix / dataframe of exogenous variables. NOT IMPLEMENTED YET
  }
  \item{holdout}{ If TRUE, the holdout sample of size h is taken from the end of the data.
  }
  \item{h}{ The length of the forecasting horizon.
  }
  \item{silent}{ If TRUE, no output or graphs are printed out.
  }
  \item{...}{ Other parameters that are passed to DON'T KNOW YET WHAT.
  }
}
\details{
The details will be here soon.
}
\value{
  The list of the following values is returned:

  \itemize{
    \item \code{persistence} - the persistence vector.
    \item \code{phi} - the value of damping parameter.
    \item \code{states} - the matrix of the components of CES. The included minimum is "level" and "potential". In the case of seasonal model the seasonal component is also included. In the case of exogenous variables the estimated coefficients for the exogenous variables are also included.
    \item \code{fitted} - the fitted values of ETS.
    \item \code{forecast} - the point forecast of ETS.
    \item \code{residuals} - the residuals of the estimated model.
    \item \code{x} - the original data.
    \item \code{ICs} - values of information criteria of the model. Includes AIC, AICc and BIC.
    \item \code{CF} - Cost function value.
    \item \code{lower} - the lower bound of prediction interval. When \code{intervals=FALSE} then NA is returned. NOT YET IMPLEMENTED
    \item \code{upper} - the higher bound of prediction interval. When \code{intervals=FALSE} then NA is returned. NOT YET IMPLEMENTED
  }
}
\references{
Hyndman, R.J., Koehler, A.B., Ord, J.K., and Snyder, R.D. (2008) Forecasting with exponential smoothing: the state space approach, Springer-Verlag. \url{http://www.exponentialsmoothing.net}.
}
\author{
Ivan Svetunkov, \email{ivan@svetunkov.ru}
}
\note{
To do:

1. Solve optimization issues.
2. Derive the prediction intervals.
3. Implement xreg.
4. Implement model selection procedure.
}

\seealso{
\code{\link[forecast]{ets}, \link[forecast]{forecast}, \link[stats]{ts}, \link[stats]{Distributions} \link[TStools]{sim.ets}}
}
\examples{
library(Mcomp)

ets2(M3$N1245$x,model="AAdN",h=12,holdout=FALSE,trace=FALSE) -> test
ets2(M3$N1245$x,model="AAdN",h=12,holdout=FALSE,trace=TRUE) -> test
ets2(M3$N2568$x,model="AMA",h=12,holdout=TRUE,trace=TRUE) -> test

}
\keyword{ exponential smoothing }
\keyword{ ETS }
\keyword{ forecasting }
\keyword{ trace likelihood }
