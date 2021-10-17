\name{review}
\alias{review}
\docType{data}
\title{review}
\description{Travel reviews}
\usage{data("review")}
\format{
  The format is:
 num [1:980, 1:10] 0.93 1.02 1.22 0.45 0.51 0.99 0.9 0.74 1.12 0.7 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:980] "User_1" "User_2" "User_3" "User_4" ...
  ..$ : chr [1:10] "Category_1" "Category_2" "Category_3" "Category_4" ...
}
\details{980 travelers' reviews of 10 different types of travel facilities in East Asia}
\source{UCI repository}
\references{Renjith, S., Sreekumar, A., and Jathavedan, M. (2018). Evaluation of partitioning clustering algorithms for processing social media data in tourism domain. 2018 IEEE Recent Advances in Intelligent Computational Systems (RAICS), 127-131.}
\examples{ 
data(review) 
## maybe str(review) ; plot(review) ... 
}
\keyword{datasets}