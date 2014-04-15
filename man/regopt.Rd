\name{regopt}
\alias{regopt}

\title{Identify regression beta using different cost functions}

\description{Identify regression beta using different cost functions.}

\usage{
regopt(y,X,constant=c(TRUE,FALSE),cost=c("MAE","MdAE","MSE","MdSE","ME","MdE"),
outplot=c(FALSE,TRUE))
}

\arguments{
  \item{y}{
    Vector of target data, can be ts object.
}
  \item{X}{
    Matrix of explanatory data, each variable is a column.
}
  \item{constant}{
    If TRUE then a constant is added to the model.
}
  \item{cost}{
    Cost function to use:
1. "MAE"  - Mean Absolute Error [Default]; 2. "MdAE" - Median Absolute Error; 3. "MSE"  - Mean Squared Error; 4. "MdSE" - Median Squared Error; 5. "ME"   - Mean Error; 6. "MdE"  - Median Error
}
  \item{outplot}{
    If TRUE plot regression fit.
}
}
\value{
\item{b}{Regression coefficients.}

}
\author{
Nikolaos Kourentzes
}
\examples{
y <- referrals
X <- c(1:length(y))
X <- as.matrix(X,nrow=length(y),ncol=1)
b <- regopt(y,X,cost="MAE",outplot=TRUE)
}