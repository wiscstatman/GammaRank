lognb <-
function( m, shape, scale )
	{
	## negative binomial mass function; log
	## make sure m >= 0, shape >0, scale > 0
	## either m may be a vector, and shape and scale should be scalers
	## or m may be a scaler, and shape and scale are vectors
	tmp1 <- lgamma(m+shape) - lgamma(shape) - lgamma(m+1)
	tmp2 <- -shape*log(scale+1)
	tmp3 <- m*( log( scale ) - log( scale + 1 ) )
	ll <- tmp1 + tmp2 + tmp3
	return(  ll )
	}

