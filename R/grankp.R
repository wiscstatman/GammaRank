grankp <-
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
		shapes <- as.integer(round(shapes))
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
	if( shapes[1] <= top )  ## corrected in version 1.1
	 {
	  Phi[1,((shapes[1]+1):(top+1)),1] <- -Inf
	 }
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
	  Beta[,j] <- tmp %*% Beta[,j+1]
	 }
	# finish
	vec <- nbrat( m=0:top, n=mhat[1],  shape=shapes[2], 
		scale=Lambda[1]/lambda[2] )  
	ok <- ( 0:top  <= shapes[1] - 1 )
	vec[!ok] <- 0
	tot <- c( Beta[,1] %*% vec )
	logp <- max( Delta[,kk-1] ) + log(tot)
	}
	}
	value <- ifelse( log.p , logp , exp(logp)  )
	return( value )
}

