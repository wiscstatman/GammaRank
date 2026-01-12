
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

#ok <- rep(FALSE, length(y))
#ok[1:25] <- TRUE
#iso.p <- grankp( post.shape[ok], post.rate[ok] )
# 25 works

## let's do a running average by size-10 window, for starters


## whole chain of length n = 162, let's go from i=6 to 157, say i +/-5
#smp <- numeric( length(y) - 10 )
#for( i in 6:157 )
# {
#  shape.sub <- y[(i-5):(i+5)]
#  smp[i] <- grankp( shapes=shape.sub, rates=rep(3,11) )
#  print(i) 
# }
### on the yearly scale there's no evidence of monotonicity!
#plot( smp )
#abline( h = -lfactorial(11), col="red" )  


### maybe if I aggregate ... say over decades...for example

## decadal data..hmm, starting at 1856...
#iceon10 <- NULL
#notdone <- TRUE
#ii <- 1:10
#d <- 1
#while( notdone )
# {
#   iceon10[d] <- sum(y[ii])
#   ii <- ii + 10
#   notdone <- !any( ii > 160 )
#   d <- d+1
# }
#smp2 <- grankp( shapes=(iceon10+1e4)[1:4], rates=rep(13, length(iceon10))[1:4]  )
## memory problem, even though aggregate data look nicer


## let's try 5 year sums

## 5-year data..hmm, starting at 1856...
iceon5 <- NULL
notdone <- TRUE
ii <- 1:5
d <- 1
while( notdone )
 {
   iceon5[d] <- sum(y[ii])
   ii <- ii + 5 
   notdone <- !any( ii > 160 )
   d <- d+1
 }
#smp3 <- grankp( shapes=(iceon5+1000), rates=rep(3, length(iceon5))  )
## to big


## Isotonic regression seems like a good testbed for gamma-rank probabilities.
## For Mendota ice-on days and a Poisson model, calculations at the annual level show very little local evidence
## for monotonicity (computing on +/- 5 year, so 11-year windows).   It might be that whatever warming process
## happens is not at the annual scale, so aggregating by simple sum, over 5 or 10 year windows retains the Poisson
## model.  The 10-year data have a pretty stong monotone signal, but computing monotone probabilities becomes
## too hard globally, for memory overflow....I try again for 5-year data, and again overflow.
## I could possibly try running windows on the 5-year or decadal data, but not sure, and have not tried yet

## Would be nice to have the code work for larger shapes, or to have an example with smaller counts
## One idea is the binomial trend test...see ../trend for more
##
