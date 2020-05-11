
# same as rl1.R, but tries a bigger shape set

source("grankp/grankp.R")
source("grankp/lognb.R")
source("grankp/nbrat.R")

n <- 20
x <- c( rnorm(n/2), rnorm(n/2, mean=0) )
groups <- c( rep(0, n/2), rep(1,n/2) )
##rr <- rank(x)

oo <- order(x, decreasing=T)  ### the anti-ranks
theta <- 2
rates <- theta^(groups[oo])
test <- grankp( shapes=as.integer(rep(3,n)), rates=rates )

## try a whole curve

theta <- seq(1/100,3, length=200 )
ll1 <- numeric( theta )
ll2 <- numeric( theta )
ll3 <- numeric( theta )
shape <- c(1,3,20)
psi <- ll1

for( i in 1:length(theta) )
 {
   rates <- (theta[i])^(groups[oo])
   ll1[i] <- grankp( shapes=as.integer( rep(shape[1],n) ), rates=rates )
   ll2[i] <- grankp( shapes=as.integer( rep(shape[2],n) ), rates=rates )
   ll3[i] <- grankp( shapes=as.integer( rep(shape[3],n) ), rates=rates )
 }
psi1 <- pbeta( theta/(1+theta), shape1=shape[1], shape2=shape[1] )
psi2 <- pbeta( theta/(1+theta), shape1=shape[2], shape2=shape[2] )
psi3 <- pbeta( theta/(1+theta), shape1=shape[3], shape2=shape[3] )

pdf( file="rl3.pdf" )
plot( psi1, ll1-max(ll1), type="l", col="black", lwd=3, ylim=c(-5,0),
		ylab="relative rank log likelihood" )
lines( psi2, ll2-max(ll2), col="red", lwd=3 )
lines( psi3, ll3-max(ll3), col="green", lwd=3 )
legend( "topright", lwd=3, col=c("black","red","green"), legend=c(
    "shape = 1", "shape = 3", "shape = 20" ) )
dev.off()

## the null case works too!
