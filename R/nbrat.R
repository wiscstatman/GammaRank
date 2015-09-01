nbrat <-
function( m, n, shape , scale )
 {
  ## returns the non-logged ratio p(m)/p(n)  of Neg binomial masses
  
  tmp1 <- lgamma(m+shape)-lgamma(n+shape)+lgamma(n+1)-lgamma(m+1)
  tmp2 <- (m-n) * ( log(scale) - log(scale+1) ) 
  ratio <- exp( tmp1 + tmp2 ) 
  ratio
 }

