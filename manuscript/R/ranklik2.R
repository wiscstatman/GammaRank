
# code to compute log rank  likelihood for location
# (no ties)
## makes plot with opposite background

library(grankp, lib.loc="~newton/Rlibs")

rl <- function( x, y, theta=seq(.01,.99,length=100), shape=5 )
 {
  ## assume no ties
  dat <- c( x, y )
  n <- length(dat)
  groups <- c( rep(0, length(x) ), rep(1, length(y) ) )
  oo <- order( dat, decreasing=TRUE )  ## anti-ranks
 
  # theta is P( X < Y ), so on the lambda scale
  b <- qbeta( theta, shape1=shape, shape2=shape )
  lambda <- (1-b)/b
  loglik <- numeric( length(theta) ) 
  for( i in 1:length(theta) )
   {
    rates <- (lambda[i])^( groups[oo] )
    loglik[i] <- grankp( shapes=as.integer(rep(shape,n)), rates=rates )
   }
  return( list( theta=theta, loglik=loglik ) ) 
 }

 ## Hollander & Wolfe (1973), 69f.
 ##    Myles Hollander and Douglas A. Wolfe (1973), _Nonparametric
  ##    Statistical Methods._ New York: John Wiley & Sons.  Pages 27-33
   ##   (one-sample), 68-75 (two-sample).
    ##  Or second edition (1999).

     ## Permeability constants of the human chorioamnion (a placental
    ##  membrane) at term (x) and between 12 to 26 weeks gestational
     ##  age (y).  The alternative of interest is greater permeability
     ##  of the human chorioamnion for the term pregnancy.

     x <- c(0.80, 0.83, 1.89, 1.04, 1.45, 1.38, 1.91, 1.64, 0.73, 1.46)
     y <- c(1.15, 0.88, 0.90, 0.74, 1.21)

gr <- as.factor( c( rep("x",length(x)), rep("y", length(y)) ) )
stripchart( c(x,y) ~ gr ) 



u <- rl(x,y, theta=seq(.05,.65,length=200) )

pdf( file="figs/pd1b.pdf"  )
yy <- u$loglik - max(u$loglik)
i <- which.min( abs( u$theta - .5 ) )
plot( u$theta, yy, type="l", lwd=3,
		xlab= ( expression(theta ==  P(X[i] < X[j])  )  ), 
		ylab="log rank likelihood - max", main="permeability data",	
		cex.axis=1.5, cex.lab=1.5 , yaxs="i" )
lines( rep( 1/2, 2 ), c( min(yy) , yy[i] ), lwd=2, col="red" )
dev.off()

## MLE at 0.28
## 2 log interval...a (0.08, 0.60)


pdf( file="figs/permb.pdf", height=3, width=6 )
plot( 0, 0, type="n", xlim=range( c(x,y) ), ylim=c(-1/2,1/2),
	axes=F , xlab="", ylab="")
abline( h=0, lwd=1 )
points( x, rep(0,length(x)), pch=19, col="red" )
points( y, rep(0,length(y)), pch=19, col="blue" )
dev.off()






