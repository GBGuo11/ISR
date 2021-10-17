\name{orange}
\alias{orange}
\docType{data}
\title{orange}
\description{orange}
\usage{data("orange")}
\format{
  The format is:
 num [1:12, 1:8] 4.79 4.58 4.71 6.58 NA ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:12] "1" "2" "3" "4" ...
  ..$ : chr [1:8] "Color.intensity" "Odor.intensity" "Attack.intensity" "Sweet" ...
}
\details{There are 19 missing values in the data set, accounting for 19.79 percent.}
\source{http://factominer.free.fr/missMDA/index.html}
\references{Josse J, Husson F (2016). missMDA: A Package for Handling Missing Values in Multivariate Data Analysis. Journal of Statistical Software, 70(1), 1–31.}
\examples{
data(orange)
## maybe str(orange) ; plot(orange) ...
}
\keyword{datasets}