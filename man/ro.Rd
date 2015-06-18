\name{ro}
\alias{ro}
\title{
Rolling Origin
}
\description{
This function makes rolling origin using the \code{data} and \code{call} passed as parameters.
}
\usage{
  ro(data,h=1,origins=1,call,value=NULL,all=FALSE)
}
\arguments{
  \item{data}{Data vector or ts object passed to the function.
}
  \item{h}{The forecasting horizon.
}
  \item{origins}{The number of rolling origins.
}
  \item{call}{The call that is passed to the function. The call must be in quotes. Example: \code{"forecast(ets(data),h)"}. Here \code{data} shows where the data and \code{h} defines where the horizon should be passed in the \code{call}.
}
  \item{value}{The variable passed by the \code{call} that contains forecasts. For example, "\code{mean}" for functions of \code{forecast} package.
}
  \item{all}{The parameter defines whether to use all the data or only that part that has the number of actuals corresponding to \code{h}.
}
}
\details{
  Function uses the defined call to produce Rolling Origin
}
\value{Function returns the following variables:
  \item{\code{forecasts} - the matrix of forecasts produced from each origin}
  \item{\code{actuals} - the matrix of actual values corresponding to the forecasts from each origin}
}
\author{
  Yves Sagaert, Ivan Svetunkov
}

\seealso{
\code{\link[TStools]{sowhat}
}
\examples{
library(forecast)
x <- rnorm(100,0,1)
call <- "forecast(ets(data),h)"

value <- "mean"
ro(x,h=5,origins=5,call,value)

value <- "$mean"
ro(x,h=5,origins=5,call,value,all=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Rolling origin }
\keyword{ Forecasting }% __ONLY ONE__ keyword per line
