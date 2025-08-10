
# 
library(rvalues)

dat <- read.csv("data2/aaSeqCDR3_138_alphaChain.csv")
dd <- cbind( dat[,3], dat[,3]+dat[,4] )
dimnames(dd) <- list( dat[,2], c("PAP","total" ) )

## from 16, 
hyps <- as.integer( c(1,299) )


B <- 10

# one order
o <- order( (dd[,1]+1/2)/(dd[,2]))
n1 <- sum(dat[,3]) ## PAP
n2 <- sum( dat[,4] ) ## ARV
N <- nrow(dat)

B <- 10  ## samples to do RB
nsim <- 100 ## number of posterior orderings
lmarg <- numeric(nsim)
bst <- numeric(B)

# repeat for posterior sample of orderings to get
# 
bigO <- matrix(NA, N, nsim)

# collect up a posterior sample of rankings
for( isim in 1:nsim )
 {
  # random draw
  lam2 <- rgamma( N, shape= hyps[1]+ dat[,4], rate=hyps[2]+n2 )
  lam1 <- rgamma( N, shape= hyps[1]+ dat[,3], rate=hyps[2]+n1 )
  bigO[,isim] <-  order( -lam1/lam2 )
 }

## for each ranking run RB, using the same seed to reduce MC

for( isim in 1:nsim )
 {
  oo <- bigO[,isim]
  set.seed(75751)  ## let's sample posterior of lambda2's  [i.e. ARV] 
  for( i in 1:B )
   {
    lam2 <- rgamma( N, shape= hyps[1]+ dat[,4], rate=hyps[2]+n2 )

    shape1 <- hyps[1] + dat[,3]

    rate1 <- lam2*(hyps[2]+n1)

    bst[i] <- grankp( shapes=shape1[oo], rates=rate1[oo], log.p=TRUE )
   }

  lmarg[isim] <- max(bst) + log( mean( exp( bst-max(bst) ) ) )
  print(isim)
 }
 
