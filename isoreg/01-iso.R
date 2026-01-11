
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
