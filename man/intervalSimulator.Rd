\name{intervalsSimulator}
\alias{intervalsSimulator}
\title{
    Simulate data from distribution and construct prediction intervals.
}
\description{
    Function simulates data, constructs prediction intervals using three methods ("Standard","2 SD" and "HM") and then returns graph and coverage, width and distance from theoretical quantiles.
}
\usage{
intervalsSimulator(randomizer=c("norm","lnorm","pois","unif","t","beta"),
                   centre=c("mean","median","ham"),
                   obs=1000,level=0.95,silent=FALSE,...)
}
\arguments{
    \item{randomizer}{ What type of random numbers generator to use.
    }
    \item{centre}{ What type of centre estimation to use (used by "HM"-based intervals only).
    }
    \item{obs}{ Number of observations to produce.
    }
    \item{level}{ Defines width of prediction interval.
    }
    \item{silent}{ If \code{TRUE}, then nothing is printed out or drawn.
    }
    \item{...}{ Parameters passed to randomizer.
    }
}
\value{
  Function returns object of class "intervalStuff", which has the following values in the list:

  \itemize{
    \item \code{coverage} - coverage for all the 3 methods used.
    \item \code{distances} - distances from upper, lower quantiles and overall.
    \item \code{width} - width of produced prediction intervals (in units).
  }
}
\author{
    Ivan Svetunkov
}
\examples{
intervalsSimulator("t","median",df=3)

intervalsSimulator("lnorm","median")
}
\keyword{ intervals }
\keyword{ prediction }
\keyword{ prediction intervals }
