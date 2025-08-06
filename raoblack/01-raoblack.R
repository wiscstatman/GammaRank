
# For Beta order probabilities, we note these are mixtures of conditionally Gamma order
# probabilities, conditioning, say, on the non-numerator Gamma variate defining the beta

# (realized Aug 6, 25 on way back from JSM Nashville)

# So analogously to Rao-Blackwellization, we estimate the order probability not be 
# averaging indicators over the joint distribution of Gamma pairs but by averaging conditiona
# Gamma order probabilities

# Let's try it

library(grankp)

# made up example


# binomial counts
# successes        5   7   10  20
# trials          20  20   40  70
