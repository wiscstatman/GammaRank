
# 
library(rvalues)

data(NBA1314)
dd <- NBA1314[,c(8,9)] ## Free throw makes vs attempts


dimnames(dd)[[1]] <- NBA1314$PLAYER

# fit via r-values

u <- rvalues(dd, family=binomial)

## hyper params

hyps <- u$aux$hypers ### alpha and beta of Beta(alpha,beta)


## idea
##
## simulate posterior directly [via independent Beta sampling]
## and use grankp to approximate the joint density per rank


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

## but I want to make simultaneous bands, possibly using grankp, which is only an approximation to
## the beta posterior

## chatgpt helps me with mean matching, given Beta(a,b), it has the same first two moments as Gamma(alpha,theta)
## where 

## shape: alpha = a^2( a+b+1)/(b*(a+b)),   rate: theta = a(a+b+1)/[b(a+b)]; gamma mean = shape/rate

## let's try computing with grankp
a <- pshapes[,1]
b <- pshapes[,2]

gshape <- a^2 * (a+b+1)/(b*(a+b) )
grate <- a*(a+b+1)/(b*(a+b))

ranklogprob <- rep(NA,B)

set.seed(75751)
for( b in 1:B )
 {
   pstar <- rbeta( nplayers, pshapes[,1], pshapes[,2] )
   o <- order(pstar,decreasing=TRUE )
   rankSim[,b] <- rank(1-pstar)  ## rank from high probabilities downward; i.e. start with best players

   # grank
   ranklogprob[b] <- grankp( shapes=gshape[o], rates=grate[o] )
print(b)

 }





