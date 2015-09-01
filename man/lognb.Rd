\name{lognb}
\Rdversion{1.1}
\alias{lognb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 utility function for grankp; log neg binom
}
\description{
 tmp
}
\usage{
lognb(m, shape, scale)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{
 count
}
  \item{shape}{ shape
} 
  \item{scale}{ scale }   
}
\details{
tmp
}
\value{
 tmp
}
\references{
tmp
}
\author{ Newton } 
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
function( m, shape, scale )
	{
	## negative binomial mass function; log
	## make sure m >= 0, shape >0, scale > 0
	## either m may be a vector, and shape and scale should be scalers
	## or m may be a scaler, and shape and scale are vectors
	tmp1 <- lgamma(m+shape) - lgamma(shape) - lgamma(m+1)
	tmp2 <- -shape*log(scale+1)
	tmp3 <- m*( log( scale ) - log( scale + 1 ) )
	ll <- tmp1 + tmp2 + tmp3
	return(  ll )
	}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
