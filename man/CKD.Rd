\name{CKD}
\alias{CKD}
\docType{data}
\title{CKD}
\description{chronic kidney disease}
\usage{data("CKD")}
\format{
  The format is:
 num [1:400, 1:18] 48 7 62 48 51 60 68 24 52 53 ...
 - attr(*, "dimnames")=List of 2
  ..$ : NULL
  ..$ : chr [1:18] "age" "bp" "sg" "al" ...
}
\details{There are 1010 missing values in the data set, accounting for 14.03 percent.}
\source{Dr.P.Soundarapandian.M.D.,D.M 
(Senior Consultant Nephrologist), 
Apollo Hospitals, 
Managiri, 
Madurai Main Road, 
Karaikudi, 
Tamilnadu, 
Indi}
\references{
Polat, H., Danaei-Mehr, H., and Cetin, A. (2017). Diagnosis of chronic kidney disease based on support vector machine by feature selection methods. Journal of Medical Systems, 41(4), 1-11.
}
\examples{
data(CKD)
## maybe str(CKD) ; plot(CKD) ...
}
\keyword{datasets}