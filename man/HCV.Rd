\name{HCV}
\alias{HCV}
\docType{data}
\title{HCV}
\description{Hepatitis C virus}
\usage{data("HCV")}
\format{
  The format is:
 num [1:615, 1:13] 1 1 1 1 1 1 1 1 1 1 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:615] "1" "2" "3" "4" ...
  ..$ : chr [1:13] "Category" "Age" "Sex" "ALB" ...
}
\details{There are 31 missing values in the data set, accounting for 0.39 percent.}
\source{UCI repository}
\references{Lichtinghagen, R., Pietsch, D., Bantel, H., Manns, M., Brand, K. and Bahr, Matthias. (2013). The Enhanced Liver Fibrosis (ELF) Score: Normal Values, Influence Factors and Proposed Cut-Off Values.. Journal of hepatology. 59. 236-242.}
\examples{
data(HCV)
## maybe str(HCV) ; plot(HCV) ...
}
\keyword{datasets}
