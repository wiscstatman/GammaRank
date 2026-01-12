
# Here's a try at computing the probability that the latent means in an regression problem are isotonic
# with app to Mendota Ice

rm(list=ls())
dat <- read.csv("MendotaIce.csv",header=TRUE)

x0 <- 2018 - dat$WINTER   # so we can consider non-decreasing
y0 <- dat$DAYS
 
ok <- !is.na(y0)   # clean data
x <- x0[ok] ; y <- y0[ok]
fit <- isoreg( x, y )

fit$x <- dat$WINTER[ok]
par( mar=c(3,3,3,0))
plot(fit, main="Ice-on days, Lake Mendota", xlab="year", ylab="days", las=1 )
legend( "bottomleft", legend="isotonic regression", lwd=1, col="red" )

### recall Gamma/Poisson conjugacy
### if lambda ~ Gamma[ alpha, beta ]; X|lambda ~ Poisson(lambda)
## then lambda | x ~ Gamma[ alpha + x, beta + 1 ]


## so shapes should be alpha + days
## rates should be beta + 1

## E(days) = E(lambda ) = alpha/beta   ### e.g. alpha = 200, beta = 2

post.shape <- as.integer( y + 200 )
post.rate <- rep(3, length(y))

## too big...let's try a small one

ok <- rep(FALSE, length(y))
ok[1:20] <- TRUE

iso.p <- grankp( post.shape[ok], post.rate[ok] )
