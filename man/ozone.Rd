\name{ozone}
\alias{ozone}
\docType{data}
\title{ozone}
\description{ozone}
\usage{data("ozone")}
\format{
  A data frame with 112 observations on the following 11 variables.
  \describe{
    \item{\code{maxO3}}{a numeric vector}
    \item{\code{T9}}{a numeric vector}
    \item{\code{T12}}{a numeric vector}
    \item{\code{T15}}{a numeric vector}
    \item{\code{Ne9}}{a numeric vector}
    \item{\code{Ne12}}{a numeric vector}
    \item{\code{Ne15}}{a numeric vector}
    \item{\code{Vx9}}{a numeric vector}
    \item{\code{Vx12}}{a numeric vector}
    \item{\code{Vx15}}{a numeric vector}
    \item{\code{maxO3v}}{a numeric vector}
  }
}
\details{There are 115 missing values in it, accounting for 9.96 percent.}
\source{http://factominer.free.fr/missMDA/index.html}
\references{Audigier, V., Husson, F., and Josse, J. (2014). A principal components method to impute missing
values for mixed data. Advances in Data Analysis and Classification, 10(1), 5-26.}
\examples{
data(ozone)
## maybe str(ozone) ; plot(ozone) ...
}
\keyword{datasets}
