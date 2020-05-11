
## First attempts at Viterbi for identifying the maximal summand

# inputs

#shapes <- c( 3, 1, 2, 2 )
#lambda <- c( 1, 2.5, 1, 1/2 )
#shapes <- rep( 3, 5 )
#shapes <- c(3,2,1)
#shapes <- c(13,11,9)
#shapes <- rep(3,10)
#shapes <- 10:1
#lambda <- 1:10
#shapes <- rep(3,100)
shapes <- c(1,10,1000)

lambda <- rep(1, length(shapes) )

## if shapes all 1 then there is a single summand!
## check for that; actually if first k-1 shapes are 1 then it's so...
## check that on input...

kk <- length( shapes )  

## make sure integer shapes (maybe with ceiling) and positive rates
## for stability it might be good to renormalize rates; e.g.
##  
lambda<- lambda/mean(lambda)   ## so they have mean 1
Lambda <- cumsum(lambda)

top <- sum( shapes[1:(kk-1)] ) - (kk-1)  ## top of the support

delta <- matrix( NA, top+1, kk-1 )
psi <- delta

lognb <- function( m, shape, scale )
	{
	## negative binomial mass function; log
	## make sure m >= 0, shape >0, scale > 0
	## m may be a vector, but shape and scale should be scalers
	tmp1 <- lgamma(m+shape) - lgamma(shape) - lgamma(m+1)
	tmp2 <- -shape*log(scale+1)
	tmp3 <- m*( log( scale ) - log( scale + 1 ) )
	ll <- tmp1 + tmp2 + tmp3
	return(  ll )
	}


Phi <- array( 0, c( top+1, top+1, kk-1 ) )
Phi[1,(1:shapes[1]),1] <- lognb( 0:(shapes[1]-1), shape=shapes[2],
			scale=Lambda[1]/lambda[2] )
Phi[1,((shapes[1]+1):(top+1)),1] <- -Inf
for( j in 2:(kk-1) )
 {
  vec <- lognb( 0:top, shape=shapes[j+1], scale=Lambda[j]/lambda[j+1] )  
  tmp2 <- matrix( rep( vec, top+1 ), top+1, top+1 , byrow=T)
  ok <- outer(0:top,0:top, function(x,y,aux=shapes[j]){ (y<= x+ aux -1 ) })
  tmp2[!ok] <- -Inf
  Phi[,, j] <- tmp2
 }

tmp1 <-  Phi[,,2] + ( Phi[1,,1] )

Delta <- matrix(NA,top+1,kk-1)  ## Do a Beta calculation if kk=2...to do..
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
ival <- mhat
mhat[kk-1] <- which.max( Delta[,kk-1] )  - 1
ival[kk-1] <- max( Delta[,kk-1] )
for( j in (kk-2):1 )
 {
  mhat[j] <- Psi[ (mhat[j+1])+1, j+1]  - 1
 }

## Next have a way to evaluate the ratio f(m_1...m_{K-1})/f( mode )
## do it one factor at a time

nbrat <- function( m, n, shape , scale )
 {
  ## returns the non-logged ratio p(m)/p(n)  of Neg binomial masses
  
  tmp1 <- lgamma(m+shape)-lgamma(n+shape)+lgamma(n+1)-lgamma(m+1)
  tmp2 <- (m-n) * ( log(scale) - log(scale+1) ) 
  ratio <- exp( tmp1 + tmp2 ) 
  ratio
 }

## Now code the backwards algorithm using the normalized summand nbrat

Beta<- matrix(1,top+1,kk-1)  ## 

for( j in (kk-2):1 )  ## make sure kk>=3
 {
  vec <- nbrat( m=0:top, n=mhat[j+1],  shape=shapes[j+2], 
		scale=Lambda[j+1]/lambda[j+2] )  
  tmp <- matrix( rep( vec, top+1 ), top+1, top+1 , byrow=T)
  ok <- outer(0:top,0:top, function(x,y,aux=shapes[j+1]){ (y<= x+ aux -1 ) })
  tmp[!ok] <- 0
  Beta[,j] <- tmp %*% Beta[,j+1]
 }
# finish
  vec <- nbrat( m=0:top, n=mhat[1],  shape=shapes[2], 
		scale=Lambda[1]/lambda[2] )  
  ok <- ( 0:top  <= shapes[1] - 1 )
  vec[!ok] <- 0
  tot <- c( Beta[,1] %*% vec )

logp <- max( Delta[,kk-1] ) + log(tot)

