# Cost weighted loss generator
lcost <- function(c) { 
	function(y,e) y*(1-c)*as.real(e<c) + (1-y)*c*as.real(e>=c)
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

# Opertator to compute the regret function for a loss
B <- function(loss) function(y,e) loss(y,e) - Lmin(loss)(y)

psi <- function(loss, c) {
	function(y,e) {
		alpha <- B(lcost(c))(y,e)
		Lmin(loss)(c) - Lmin
	}
}

# Example use of expressions and automatic differentiation 
eval(deriv(expression(3*x^2 + 2), "x"), environment(x <- 0.3))[1]