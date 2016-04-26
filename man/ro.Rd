\name{ro}
\alias{ro}
\title{
Rolling Origin
}
\description{
This function makes rolling origin using the \code{data} and \code{call} passed as parameters.
Currently function can do all of that in parallel, but it needs \code{foreach} and either \code{doMC} (Linux only), \code{doParallel} or \code{doSNOW} packages installed.
}
\usage{
  ro(data,h=10,origins=10,call,value=NULL,
     ci=FALSE,co=FALSE,silent=FALSE,parallel=FALSE)
}
\arguments{
  \item{data}{Data vector or ts object passed to the function.
  }
  \item{h}{The forecasting horizon.
  }
  \item{origins}{The number of rolling origins.
  }
  \item{call}{The call that is passed to the function. The call must be in quotes. Example: \code{"forecast(ets(data),h)"}. Here \code{data} shows where the data and \code{h} defines where the horizon should be passed in the \code{call}. Some hidden parameters can also be specified in the call. For example, parameters \code{counti}, \code{counto} and \code{countfull} are used in the inner loop and can be used for the regulation of exogenous variables sizes. See examples for the appropriate command.
  }
  \item{value}{The variable passed by the \code{call} that contains forecasts. For example, "\code{mean}" for functions of \code{forecast} package. Can also be a vector of variables. See examples for the details. If the parameter is \code{NULL}, then all the values from the call are returned (could be really messy!).
  }
  \item{ci}{The parameter defines if the in-sample window size should be constant.
  }
  \item{co}{The parameter defines whether the holdout sample window size should be constant. If \code{TRUE} the rolling origin will stop when less than \code{h} observations are left.
  }
  \item{silent}{If TRUE, nothing is printed out in the console.
  }
  \item{parallel}{If TRUE, then the rolling origin is done in parallel. WARNING! Packages \code{foreach} and either \code{doMC} (Linux only), \code{doParallel} or \code{doSNOW} are needed in order to run the function in parallel.
  }
}
\details{
  Function uses the defined call to produce Rolling Origin
}
\value{Function returns the following variables:
\itemize{
  \item{\code{actuals} - the matrix of actual values corresponding to the produced forecasts from each origin.}
  \item{value - the matrices with the produced data from each origin. Name of each matrix corresponds to the names in the parameter \code{value}.}
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
ourcall <- "forecast(ets(data),h=h)"

ourvalue <- "mean"
ro(x,h=5,origins=5,ourcall,ourvalue)

ourvalue <- "$mean"
ro(x,h=5,origins=5,ourcall,ourvalue,ci=TRUE)

ourvalue <- c("$mean","lower[,1]","upper[,1]","lower[,2]","upper[,2]")
ro(x,h=5,origins=20,ourcall,ourvalue,ci=TRUE,co=TRUE)

#### The following code will give exactly the same result as above,
#### but computed in parallel using all but 1 core of CPU:
# ro(x,h=5,origins=20,ourcall,ourvalue,ci=TRUE,co=TRUE,parallel=TRUE)

#### A more complicated example using the for loop and
#### first 5 series from M3-competition data

library(Mcomp)
## Make a subset from M3 data
M3.yearly <- subset(M3,1)

## Form an array for the forecasts we will produce
M3.forecasts <- array(NA,c(6,4,5))

## Define models that need to be used for each series
ourmodels <- c("AAN","ANN","MNN","MMN","MNN")

## This call uses specific models for each time series
ourcall <- "forecast(ets(data),model='ourmodels[i]',h=h)"
ourvalue <- "mean"

## Start the loop. The important thing here is to use the same variable i as in ourcall.
for(i in 1:5){
    ourdata <- ts(c(M3.yearly[[i]]$x,M3.yearly[[i]]$xx))
    M3.forecasts[,,i] <- ro(data=ourdata,h=6,origins=4,call=ourcall,
                            value=ourvalue,co=TRUE,silent=TRUE)$mean
}

## M3.forecasts array now contains rolling origin forecasts from specific
## models for first 5 time series in M3 yearly data.
M3.forecasts


##### An example with exogenous variables
x <- rnorm(100,0,1)
xreg <- rnorm(100,0,1)

## counti is used to define in-sample size of xreg, counto - the size of holdout sample
ourcall <- "forecast(auto.arima(data,xreg=xreg[counti]),h=h,xreg=xreg[counto])"
ourvalue <- "mean"

ro(x,h=5,origins=5,ourcall,ourvalue)

## countfull is used to take xreg of the size corresponding to all sample on each iteration
ourcall <- "es(data,'ANN',xreg=xreg[countfull],h=h,silent=TRUE)"
ourvalue <- "forecast"

ro(x,h=5,origins=5,ourcall,ourvalue)
}
\keyword{ Rolling origin }
\keyword{ Forecasting }