\name{nbrat}
\Rdversion{1.1}
\alias{nbrat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 a utility function for grankp.. ratio of negative binomials
}
\description{
 ratio of nb
}
\usage{
nbrat(m, n, shape, scale)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{
 tmp
}
  \item{n}{
 tmp
}
  \item{shape}{
 tmp
}
  \item{scale}{
 tmp
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ tmp }

 
\references{
%% ~put references to the literature/web site here ~
}
\author{
M A Newton
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function( m, n, shape , scale )
 {
  ## returns the non-logged ratio p(m)/p(n)  of Neg binomial masses
  
  tmp1 <- lgamma(m+shape)-lgamma(n+shape)+lgamma(n+1)-lgamma(m+1)
  tmp2 <- (m-n) * ( log(scale) - log(scale+1) ) 
  ratio <- exp( tmp1 + tmp2 ) 
  ratio
 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
