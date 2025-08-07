
# 
library(rvalues)

data(NBA1314)
dd <- NBA1314[,c(8,9)] ## Free throw makes vs attempts


dimnames(dd)[[1]] <- NBA1314$PLAYER

# fit via r-values

u <- rvalues(dd, family=binomial)

## hyper params

hyps <- u$aux$hypers ### alpha and beta of Beta(alpha,beta)
hyps <- as.integer( hyps )  ## for grankp


## idea (an improvement over what I started in `../HPD/examp1.R` 
##
## simulate posterior directly [via independent Beta sampling]
## and use grankp to approximate the joint density per rank ... BUT HERE USE RAO-BLACKWELLIZED rank probs
## the approximation I tried in ../HPD/examp1.R doesn't work for large shapes...back to square 1


B <- 10
stats <- u$main ## rvalue output, r-value sorted
nplayers <- nrow(stats)

rankSim <- matrix(NA, nplayers, B )  ## rows per player, columns per posterior sample
pshapes <- cbind( stats$xx + hyps[1], stats$nn-stats$xx + hyps[2] )

set.seed(75751)
for( b in 1:B )
 {
   pstar <- rbeta( nplayers, pshapes[,1], pshapes[,2] )
   rankSim[,b] <- rank(1-pstar)  ## rank from high probabilities downward; i.e. start with best players
 }

## OK, so the above gives me a posterior sample of the ranks

## but I want to make simultaneous bands, possibly using grankp, 

# let's try grankp RB just to evaluate log.p of the top ranking

## 

#nr <- 100
#bst <- rep(NA,nr)
#kp <- rep(FALSE, nplayers)
#kp[1:10] <- TRUE
###kp[1:4] <- TRUE
#for( i in 1:nr )
# {
#  u <- rgamma( nplayers[kp], shape=pshapes[kp,2])
#  bst[i] <- grankp( shapes=pshapes[kp,1], rates=u, log.p=TRUE )
#  print(i)
# }
#
#lmarg <- max(bst) + log( mean( exp( bst-max(bst) ) ) )

##
### pass 1: sum(kp)=25...top 25
## seems ok; lmarg = -74.1; while max(bst) = -69.5, so we surely get something by averaging
## that's the MC estimate of the marginal mass of the ordering of the r-value ranked top 25
## 
## but it's weird...the log( 1/25! ) = -58, which is quite a bit higher than -74...hmm
## maybe test it with top 10 instead, for speed

### pass 10; sum(kp)=10...top 10...; now lmarg = -18; max(bst) = -14.2, but log(10!) = 15.1...so lmarg is still oddly small

### pass 4; sum(kp)=4
## this one ok
#> max(bst)
#[1] -0.06520867
#> lmarg
#[1] -2.526301
#> abline( v=lmarg )
#> hist( bst, 30 )
#> abline( v=lmarg )
#> exp(lmarg)
#[1] 0.07995426
#> 1/factorial(4)
#[1] 0.04166667


## I should confirm in some basic cases that the total mass over all k! perms equals 1, say where mass is exp(lmarg)
#  computed by RB

library(partitions)

##pp <- partitions::perms(4)  ## the 24 perms
pp <- partitions::perms(7)  ## the go big


lmarg <- numeric(24)
nr <- 100
bst <- rep(NA,nr)
kp <- rep(FALSE, nplayers)
#kp[1:4] <- TRUE
kp[1:7] <- TRUE
#
nplay <- nplayers[kp]
s1 <- pshapes[kp,1]
s2 <- pshapes[kp,2]
for( j in 1:ncol(pp) )
 {
  set.seed(75751)  ## let's use the same Gamma's
  rr <- pp[,j]

  for( i in 1:nr )
   {
    u <- rgamma( nplay, shape=s2 )
    bst[i] <- grankp( shapes=s1[rr], rates=u[rr], log.p=TRUE )
   }
  lmarg[j] <- max(bst) + log( mean( exp( bst-max(bst) ) ) )
  print(j)
 }
## sum( exp(lmarg) ) = 1 ... it works, well at least with k=4!



