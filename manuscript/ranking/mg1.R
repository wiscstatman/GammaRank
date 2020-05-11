
# Mike G asks about confidence in rankings. (Jan 2011)
# e.g., he has sample of n=20 types from a population, and sees
# 10 A's , 5 B's , 2 C's and one each of D,E, F, and no G's orH's.

# 

library( grankp, lib.loc="~newton/Rlibs" )
xx <- as.integer( c(10,5,2,1,1,1,0,0) + 1 ) 
pp <- 6*2*grankp( xx , log.p=FALSE )
