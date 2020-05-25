
library(grankp, lib.loc="~newton/Rlibs")


## code to look at parameters where modal ordering differs from posterior mean ordering

# example (found by Nick when we were revising the rvalue paper)

t <- c(0.2, 1.9, 0.8)
mu <- c(0.85, 1.55, 0.95)/( 1 + t) 
t <- c(0.2, 1.9, 0.8)
sigma2 <- t/(1+t)

n <- 2000

ss <- as.integer( round(n/sigma2) )
rr <- ss*(1-mu/sqrt(n)) 

test <- numeric(6)
ord <- c(1,2,3)
test[1] <- grankp( ss[ord], rr[ord], log.p=FALSE )
names(test)[1] <- paste( ord , collapse="-" )

ord <- c(1,3,2)
test[2] <- grankp( ss[ord], rr[ord], log.p=FALSE )
names(test)[2] <- paste( ord , collapse="-" )

ord <- c(2,1,3)
test[3] <- grankp( ss[ord], rr[ord], log.p=FALSE )
names(test)[3] <- paste( ord , collapse="-" )

ord <- c(2,3,1)
test[4] <- grankp( ss[ord], rr[ord], log.p=FALSE )
names(test)[4] <- paste( ord , collapse="-" )

ord <- c(3,1,2)
test[5] <- grankp( ss[ord], rr[ord], log.p=FALSE )
names(test)[5] <- paste( ord , collapse="-" )

ord <- c(3,2,1)
test[6] <- grankp( ss[ord], rr[ord], log.p=FALSE )
names(test)[6] <- paste( ord , collapse="-" )

## simulation check
B <- 10^7
x <- rnorm( B, mean=mu[1], sd=sqrt(sigma2[1]) )
y <- rnorm( B, mean=mu[2], sd=sqrt(sigma2[2]) )
z <- rnorm( B, mean=mu[3], sd=sqrt(sigma2[3]) )

sim <- test
sim[1] <- mean( (x>y) & (y>z) )
sim[2] <- mean( (x>z) & (z>y) )
sim[3] <- mean( (y>x) & (x>z) )
sim[4] <- mean( (y>z) & (z>x) )
sim[5] <- mean( (z>x) & (x>y) )
sim[6] <- mean( (z>y) & (y>x) )
