\name{PM2.5}
\alias{PM2.5}
\docType{data}
\title{PM2.5}
\description{Beijing PM2.5}
\usage{data("PM2.5")}
\format{
  The format is:
 num [1:43824, 1:12] 2010 2010 2010 2010 2010 2010 2010 2010 2010 2010 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:43824] "1" "2" "3" "4" ...
  ..$ : chr [1:12] "year" "month" "day" "hour" ...
}
\details{It records 43824 daily measurements on 12 variables and there are 2067 missing values on 2067 measurements, accounting for 0.00393.}
\source{UCI repository}
\references{X. Liang, T. Zou, B. Guo, S. Li, H. Zhang, S. Zhang, H. Huang, and S. Chen. Assessing Beijing's PM2.5 pollution: severity, weather impact, APEC and winter heating. Proceedings of the Royal Society A, 471(2182):1–20, 2015.}
\examples{
data(PM2.5)
## maybe str(PM2.5) ; plot(PM2.5) ...
}
\keyword{datasets}
