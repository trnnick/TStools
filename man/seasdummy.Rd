\name{seasdummy}
\alias{seasdummy}

\title{Create seasonal dummy variables.}

\description{Create binary or trigonometric seasonal dummies.}

\usage{
seasdummy(n,m=NULL,y=NULL,type=c("bin","trg"),full=c(FALSE,TRUE))
}

\arguments{
  \item{n}{
    Number of observations to create.
    }
  \item{m}{
    Seasonal periodicity. If NULL it will take the information from the provided time series (y argument). 
    }
   \item{y}{
    This is an optional time series input that can be used to get seasonal periodicity (m) and the start point. 
    }
  \item{type}{
    Type of seasonal dummies to create. Can be "bin" for binary and "trg" for trigonometric.
    }
  \item{full}{
    If full is TRUE then keeps the m-th dummy that is co-linear to the rest.
    }
}
\value{
\item{x}{Array with seasonal dummies.}
}
\author{
Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}
}
\examples{
x <- seasdummy(24,12)
}