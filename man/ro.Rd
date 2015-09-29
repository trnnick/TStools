\name{ro}
\alias{ro}
\title{
Rolling Origin
}
\description{
This function makes rolling origin using the \code{data} and \code{call} passed as parameters.
}
\usage{
  ro(data,h=1,origins=1,call,value=NULL,
    ci=FALSE,co=FALSE,silent=FALSE)
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
  \item{value}{The variable passed by the \code{call} that contains forecasts. For example, "\code{mean}" for functions of \code{forecast} package. Can also be a vector of variables. See examples for the details.
}
  \item{ci}{The parameter defines if the in-sample window size should be constant.
}
  \item{co}{The parameter defines whether the holdout sample window size should be constant. If \code{TRUE} the rolling origin will stop when less than \code{h} observations are left.
}
  \item{silent}{If TRUE, no output is given.
}
}
\details{
  Function uses the defined call to produce Rolling Origin
}
\value{Function returns the following variables:
\itemize{
  \item{\code{actuals} - the matrix of actual values corresponding to the produced forecasts from each origin.}
  \item{value - the matrices with the produced data from each origin. Name of each matrix corresponds to the names in value.}
  }
}
\author{
  Yves Sagaert, Ivan Svetunkov
}

\seealso{
\code{\link[TStools]{sowhat}}
}
\examples{
library(forecast)
x <- rnorm(100,0,1)
call <- "forecast(ets(data),h=h)"

value <- "mean"
ro(x,h=5,origins=5,call,value)

value <- "$mean"
ro(x,h=5,origins=5,call,value,ci=TRUE)

value <- c("$mean","lower[,1]","upper[,1]","lower[,2]","upper[,2]")
ro(x,h=5,origins=5,call,value,ci=TRUE,co=TRUE)
}
\keyword{ Rolling origin }
\keyword{ Forecasting }