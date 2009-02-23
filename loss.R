setwd('~/code/surrogate/')
source("Deriv.R")

# Cost weighted loss generator
lcost <- function(c) {
	function(y,e) ifelse(e<c, y*(1-c), (1-y)*c)
}

# 0-1 Loss
l01 <- function(y,e) 2*lcost(0.5)(y,e)

# Square loss
lsquare <- function(y,e) y*(1-e)^2 + (1-y)*e^2

#Log loss
llog <- function(y,e) -y*log(e) - (1-y)*log(1-e)

# Exponential Loss
lexp <- function(y,e) y*sqrt((1-e)/e) + (1-y)*sqrt(e/(1-e)) 

# Operator that computes point-wise Bayes loss from loss
Lmin <- function(loss) function(y) loss(y,y)

# Operator to compute the regret function for a loss
B <- function(loss) function(y,e) loss(y,e) - Lmin(loss)(y)

psi <- function(loss, c, alpha) {
	minloss <- Lmin(loss)
	dminloss <- Deriv.function(minloss)
	minloss(c) - minloss(c-alpha) + alpha*dminloss(c)
#	if((c-alpha < 0) || (c-alpha > 1)) {
#		return(0)
#	} else {
#		return(minloss(c) - minloss(c-alpha) + alpha*dminloss(c))
#	}
}

bound <- function(loss, c) {
	Vectorize(function(alpha) {
		min(psi(loss, c, alpha), psi(loss, c, -alpha))
	})
}
