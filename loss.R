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
B <- function(loss) Vectorize(function(y,e) loss(y,e) - Lmin(loss)(y))

psi <- function(loss, c, alpha) {
	minloss <- Lmin(loss)
	dminloss <- Deriv.function(minloss)
	minloss(c) - minloss(c-alpha) + alpha*dminloss(c)
}

bound <- function(loss, c) {
	function(y,e) {
		alpha <- B(lcost(c))(y,e)
		max(psi(loss, c, alpha), psi(loss, c, -alpha))
	}
}

plotBounds <- function(loss, c, est, title=deparse(body(loss))) {
	Bc <- B(lcost(c))
	Bloss <- B(loss)
	plot(function(y) Bloss(y,est), 0, 1, ylim=c(-0.5,1), lty=1, main=title, ylab='Regret', xlab='η')
	plot(function(y) psi(loss, c, Bc(y,est)), add=TRUE, col='red')
	plot(function(y) psi(loss, c, -Bc(y,est)), add=TRUE, col='red', lty=2)
}

invBound <- function(loss, c) {
	function(y,e) {
		b <- B(loss)(y,e)
		f <- function(alpha) max(psi(loss, c, alpha), psi(loss, c, -alpha)) - b
		uniroot(f, c(0,min(c,1-c)-0.001), maxiter=100)$root
	}
} 

# The following code tests whether estimates other than the mean make sense
# mu is a Chernoff bound for the true probability t given an estimate e built from
# n 0-1 samples drawn from Binomial(t, 1-t)
mu <- function(n,e) function(t) 2*exp(-2*n*abs(t-e)^2)

# Given a loss l, the number of sample n and the empirical estimate e,
# this returns a function that, for a new estimate s, computes an upper bound 
# on the expected true risk for l under the Chernoff bound given by mu
risk <- function(l,n,e) function(s) integrate(function(t) l(t,s)*(mu(n,e)(t)), 0, 1)$value

# Plot the bound on the expected true risk given the estimate 0.2 and the risk assuming the 
# estimate is true
plot(Vectorize(risk(lexp,5,0.2)), 0.1, 0.6, ylim=c(0.5,1))
plot(Vectorize(function(x) lexp(0.2,x)), 0.1, 0.6, add=TRUE, col=2)

# Compute the estimate that minimises the upper bound on the true risk for exponential
# loss when there are 10 samples with an empirical estimate of 0.2
optimize(risk(lexp,10,0.2), c(0.1, 0.9))

n = 5
e = 0.25
optimize(risk(lexp,n,e), c(0.01, 0.99))
optimize(risk(lsquare,n,e), c(0.01, 0.99))
optimize(risk(llog,n,e), c(0.01, 0.99))

regret <- function(l,n,e) function(s) integrate(function(t) ( l(t,s) - Lmin(l)(t) )*(mu(n,e)(t)), 0, 1)$value

optimize(regret(lexp,n,e), c(0.01, 0.99))
optimize(regret(lsquare,n,e), c(0.01, 0.99))
optimize(regret(llog,n,e), c(0.01, 0.99))

plot(Vectorize(regret(lexp,n,e)), 0.1, 0.6, ylim=c(0,0.5))
plot(Vectorize(regret(lsquare,n,e)), 0.1, 0.6, add=TRUE, col=2)
plot(Vectorize(regret(llog,n,e)), 0.1, 0.6, add=TRUE, col=3)
plot(Vectorize(regret(l01,n,e)), 0.1, 0.6, add=TRUE, col=4)


# Some plots of L(e,e') on the unit square
x <- seq(0.1,0.9,len=30)
y <- x
z <- outer(x, y, llog)
persp(x,y,z, theta=0, phi=0) -> res
lines(trans3d(x,x,z=llog(x,x), pm=res), col=3, lwd=3)
