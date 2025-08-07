
# For Beta order probabilities, we note these are mixtures of conditionally Gamma order
# probabilities, conditioning, say, on the non-numerator Gamma variate defining the beta

# (realized Aug 6, 25 on way back from JSM Nashville)

# So analogously to Rao-Blackwellization, we estimate the order probability not be 
# averaging indicators over the joint distribution of Gamma pairs but by averaging conditiona
# Gamma order probabilities

# Let's try it

library(GammaRank)

# made up example


# binomial counts
#                  A   B    C
# successes        5  10   20
# trials          10  20   40 
# say with a Beta(1,1) prior so
#
# the posteriors are Beta_A(6,12); Beta_B(11,22); and Beta_C(21,42)
# there are 6 orders...I'm curious what the order probs are

e1 <- function(nsim=1000,AA=c(6,6),BB=c(11,11),CC=c(21,21) )
 {
	pord <- matrix(NA,nsim,6)
	dimnames(pord)[[2]] <- c("ABC","ACB", "BAC", "CAB", "BCA", "CBA" )
	# let's consider the order as decreasging
	## e.g. P[ABC] = P[ A > B > C ] = E{ P[ X_A > X_B > X_C | 2nd gammas u ] }
	## where X_C = Gamma[21,1]/u[3]  u[3] = realization of Gamma[21,1]

	A1 <- AA[1]; A2 <- AA[2]
	B1 <- BB[1]; B2 <- BB[2]
	C1 <- CC[1]; C2 <- CC[2]

	for( i in 1:nsim )
	 { 
	  # sample the 3 second Gammas
	  u <- rgamma( 3, shape=c(A2,B2,C2) )  ## A , B, C
	  # conditionally on U, the Beta Ratios are 
	  # Gamma( A1, u[1] ), Gamma(B1, u[2]), Gamma(C1 ,u[3] )	
  
	  pord[i,1] <- grankp( shapes=c(A1,B1,C1), rates=u , log.p=FALSE)	
	  pord[i,2] <- grankp( shapes=c(A1,C1,B1), rates=u[c(1,3,2)], log.p=FALSE )	
	  pord[i,3] <- grankp( shapes=c(B1,A1,C1), rates=u[c(2,1,3)], log.p=FALSE )	
	  pord[i,4] <- grankp( shapes=c(C1,A1,B1), rates=u[c(3,1,2)], log.p=FALSE )	
	  pord[i,5] <- grankp( shapes=c(B1,C1,A1), rates=u[c(2,3,1)], log.p=FALSE )	
	  pord[i,6] <- grankp( shapes=c(C1,B1,A1), rates=u[c(3,2,1)], log.p=FALSE )	
	 }

	phat.raoblack <- colMeans(pord)
	return(phat.raoblack)
	}

## let's do it by simple MC
e2 <- function(nsim=1000,AA=c(6,6),BB=c(11,11),CC=c(21,21) )
 {
	Iord <- matrix(0,nsim,6)
	# let's consider the order as decreasging
	## e.g. P[ABC] = P[ A > B > C ] = E{ P[ X_A > X_B > X_C | 2nd gammas u ] }
	## where X_C = Gamma[21,1]/u[3]  u[3] = realization of Gamma[21,1]

	A1 <- AA[1]; A2 <- AA[2]
	B1 <- BB[1]; B2 <- BB[2]
	C1 <- CC[1]; C2 <- CC[2]

	for( i in 1:nsim )
	 { 
	  # sample the 3 second Gammas
  	  u <- rgamma( 3, shape=c(A2,B2,C2) )  ## A , B, C
	  w <- rgamma( 3, shape=c(A1,B1,C1) ) 
	  rr <- w/u
	  if( rr[1] > rr[2] & rr[2] > rr[3] ){ Iord[i,1] <- 1  } 
	  if( rr[1] > rr[3] & rr[3] > rr[2] ){ Iord[i,2] <- 1  } 
	  if( rr[2] > rr[1] & rr[1] > rr[3] ){ Iord[i,3] <- 1  } 
	  if( rr[3] > rr[1] & rr[1] > rr[2] ){ Iord[i,4] <- 1  } 
	  if( rr[2] > rr[3] & rr[3] > rr[1] ){ Iord[i,5] <- 1  } 
	  if( rr[3] > rr[2] & rr[2] > rr[1] ){ Iord[i,6] <- 1  } 
 	}
	phat.simple <- colMeans(Iord)
	return(phat.simple)
	}


## Now repeat both to check MC error variation
nb <- 20
phat.simple <- matrix( NA, nb, 6 )
dimnames(phat.simple)[[2]] <- c("ABC","ACB", "BAC", "CAB", "BCA", "CBA" )
phat.rb <- phat.simple

for( b in 1:20 )
 {
  phat.rb[b,] <- e1()
  phat.simple[b,] <- e2()
  print(b)
 }

apply( phat.rb, 1, sd )

apply( phat.simple,1 sd )  ## check the means first, and they're on target
##  Rao-blackwellization  works!  RB has about 1/2 the sd for the same number draws nsim!
## ...maybe other cases the sd improvement will be even more


