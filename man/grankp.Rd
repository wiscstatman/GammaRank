\name{grankp}
\Rdversion{1.1}
\alias{grankp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to compute the rank probabilities for independent but not
identically distributed Gamma random variables.
}
\description{
With \eqn{Z_1, Z_2, \ldots , Z_k} equal to independent Gamma-distributed
random variables with positive integer shapes \eqn{a_1, a_2, \ldots, a_k}
 and positive rates \eqn{\lambda_1, \lambda_2, \ldots, \lambda_k}, 
 returned is the (log) probability 
 \eqn{P[Z_1> Z_2 > \ldots, > Z_k ]}
}
\usage{
grankp(shapes, rates = rep(1, length(shapes)), log.p = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{shapes}{
 A vector of positive integer shape parameters.
}
  \item{rates}{
 A vector of positive rate parameters
}
  \item{log.p}{
 logical indicating whether or not to return natural logarithm
}
}
\details{
When $k=2$ the computation is a simple Beta cumulative probability;
when $k>2$ the code uses a negative-binomial representation of the
event (Newton and Chung, 2010, Annals of Statistics, 38, 3217-3244), and a hidden-markov-chain-style
Viterbi and backward recursion (Newton, in prep).  The calculation
is central to clustering by gamma-ranking. 
}
\value{
log probability  of the ordering event
}
\references{ 
 \url{http://projecteuclid.org/euclid.aos/1284988405}
\author{
Michael Newton \url{http://www.stat.wisc.edu/~newton}
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ }
\examples{


n <- 7
grankp( as.integer( rep(5,n) ) )

# compare to
log( 1/factorial(n) )

# In spring 2008, a random sample of n=30 votors were each asked 
# if they prefer Obama, McCain, or H Clinton, and the respective
# favorite counts were 12, 10, and 8.  Under a flat Dirichlet
# prior for the population proportions (p_O, p_M, p_C), the
# posterior probability that the empirical ordering is the population 
# ordering is

grankp( as.integer( c(13,11,9) ), log.p=FALSE )

## The function is currently defined as
function( shapes, rates=rep(1,length(shapes)), log.p=TRUE )
 {
	# computes log Prob[ Z_1 > ... > Z_kk ]
        # where Z_i ~ Gamma[ shape[i], rate[i] 

	kk <- length( shapes )  
	lambda <- rates  ## renamed for consistency with paper

	## check input
	if( length(rates) != length(shapes) )
		{ 
		stop( "length(rates) != length(shapes)" )
		}
	if( any( rates <= 0 ) )
		{
		stop( "All rates must be non-negative" )
		}
	if( any( is.na(shapes) ) | any( is.na(rates) ) )
		{
		stop("missing values not allowed")
		}
	if( !is.integer(shapes) )
		{
		warning( "shapes coerced to integers as possible" )
		shapes <- as.integer(shapes)
		}
	if( any( shapes ) <= 0 )
		{
		stop( "shapes must be positive integers" )
		}
	lambda<- lambda/mean(lambda)   ## so they have mean 1
	Lambda <- cumsum(lambda)  ## cumulative values
	
	if( kk == 2 )
	{
	## do a Beta calculation
	crit <- lambda[2]/Lambda[2]
	logp <- pbeta( crit, shape1=shapes[1], shape2=shapes[2], 
		lower.tail=FALSE, log.p=TRUE )
	}
	if( kk >= 3 )
	{
	simple <- all( shapes[1:(kk-1)] == 1 )
	if( simple )
	{
	## there is a single summand
	logp <-  sum( lognb(0, shape=shapes[2:kk], 
			scale=(Lambda[1:(kk-1)]/lambda[2:kk]) ) )
		## note, the last shape, shapes[kk], might exceed 1...
	}
	if( !simple )
	{
	## do the full HMM-style calculation
	top <- sum( shapes[1:(kk-1)] ) - (kk-1)  ## top of the support

	## Viterbi to find modal summand
	delta <- matrix( NA, top+1, kk-1 )
	psi <- delta

	Phi <- array( 0, c( top+1, top+1, kk-1 ) )
	Phi[1,(1:shapes[1]),1] <- lognb( 0:(shapes[1]-1), shape=shapes[2],
			scale=Lambda[1]/lambda[2] )
	Phi[1,((shapes[1]+1):(top+1)),1] <- -Inf
	for( j in 2:(kk-1) )
	 {
	  vec <- lognb( 0:top, shape=shapes[j+1], scale=Lambda[j]/lambda[j+1] ) 
	tmp2 <- matrix( rep( vec, top+1 ), top+1, top+1 , byrow=TRUE)
	ok <- outer(0:top,0:top, function(x,y,aux=shapes[j]){ (y<= x+ aux -1)})
	tmp2[!ok] <- -Inf
	Phi[,, j] <- tmp2
	}
	tmp1 <-  Phi[,,2] + ( Phi[1,,1] )

	Delta <- matrix(NA,top+1,kk-1) 
	Psi<- matrix(NA,top+1,kk-1)
	Delta[,2] <- apply( tmp1, 2, max )  ## ok if kk=3
	Psi[,2] <- apply( tmp1, 2, which.max )
	if( kk >= 4 )
	 {
	  for( j in 3:(kk-1) )
	   {
	   tmp3 <- matrix(rep( Delta[,j-1], top+1 ), top+1, top+1 ) + Phi[,,j] 
	   Delta[,j] <- apply( tmp3, 2, max )
	   Psi[,j] <- apply( tmp3, 2, which.max )
	  }
	 }

	## backtrack 
	mhat <- numeric( kk - 1)
	mhat[kk-1] <- which.max( Delta[,kk-1] )  - 1
	for( j in (kk-2):1 )
	 {
	  mhat[j] <- Psi[ (mhat[j+1])+1, j+1]  - 1
	 }
	## end Viterbi step

	## Backwards algorithm using the normalized summands
	Beta<- matrix(1,top+1,kk-1)  ## 
	for( j in (kk-2):1 )  ## make sure kk>=3
	 {
	  vec <- nbrat( m=0:top, n=mhat[j+1],  shape=shapes[j+2], 
		scale=Lambda[j+1]/lambda[j+2] )  
	  tmp <- matrix( rep( vec, top+1 ), top+1, top+1 , byrow=TRUE)
	  ok <- outer(0:top,0:top, 
		function(x,y,aux=shapes[j+1]){ (y<= x+ aux -1 ) })
	  tmp[!ok] <- 0
	  Beta[,j] <- tmp \%*\% Beta[,j+1]
	 }
	# finish
	vec <- nbrat( m=0:top, n=mhat[1],  shape=shapes[2], 
		scale=Lambda[1]/lambda[2] )  
	ok <- ( 0:top  <= shapes[1] - 1 )
	vec[!ok] <- 0
	tot <- c( Beta[,1] \%*\% vec )
	logp <- max( Delta[,kk-1] ) + log(tot)
	}
	}
	value <- ifelse( log.p , logp , exp(logp)  )
	return( value )
  }
}
\keyword{distribution}
