
# inputs

#shapes <- rep( 6, 5 )
#lambda <- rep(1, 5 )
#shapes <- c(13,11,9)
#shapes <- c(3,2,1)
#shapes <- c(9,11,13)
#shapes <- c(100,100,100)
#shapes <- rep(3, 10)
#shapes <- 10:1
lambda <- 1:10

#shapes <- rep(3,100)

shapes <- c(1,10,1000)
lambda <- rep(1,length(shapes))


kk <- length( shapes )  


## make sure integer shapes (maybe with ceiling) and positive rates
## for stability it might be good to renormalize rates; e.g.
##  
lambda<- lambda/mean(lambda)   ## so they have mean 1
Lambda <- cumsum(lambda)
top <- sum( shapes[1:(kk-1)] ) - (kk-1)  ## top of the support


nb <- function( m, shape, scale )
	{
	## negative binomial mass function; log
	## make sure m >= 0, shape >0, scale > 0
	## m may be a vector, but shape and scale should be scalers
	tmp1 <- lgamma(m+shape) - lgamma(shape) - lgamma(m+1)
	tmp2 <- -shape*log(scale+1)
	tmp3 <- m*( log( scale ) - log( scale + 1 ) )
	ll <- tmp1 + tmp2 + tmp3
	return( exp( ll ) )
	}

## Now code the backwards algorithm using the normalized summand nbrat

Beta<- matrix(1,top+1,kk-1)  ## 

for( j in (kk-2):1 )  ## make sure kk>=3
 {
  vec <- nb( m=0:top,   shape=shapes[j+2], 
		scale=Lambda[j+1]/lambda[j+2] )  
  tmp <- matrix( rep( vec, top+1 ), top+1, top+1 , byrow=T)
  ok <- outer(0:top,0:top, function(x,y,aux=shapes[j+1]){ (y<= x+ aux -1 ) })
  tmp[!ok] <- 0
  Beta[,j] <- tmp %*% Beta[,j+1]
 }
# finish
  vec <- nb( m=0:top,   shape=shapes[2], 
		scale=Lambda[1]/lambda[2] )  
  ok <- ( 0:top  <= shapes[1] - 1 )
  vec[!ok] <- 0
  tot <- c( Beta[,1] %*% vec )

logp <- log(tot)
## ok, but doesn't use the viterbi normalization yet; seems to work
