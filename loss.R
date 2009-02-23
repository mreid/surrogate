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
	Vectorize(function(alpha) {
		min(psi(loss, c, alpha), psi(loss, c, -alpha))
	})
}

bound2 <- function(loss, c) {
	Bc <- B(lcost(c))
	function(y,e) {
		alpha = Bc(y,e)
		max(psi(loss, c, alpha), psi(loss,c,-alpha))
	}
}

# Plot some bounds as a function of alpha
bsq <- bound(lsquare, 0.5)
bexp <- bound(lexp, 0.5)
blog <- bound(llog, 0.5)

plot(bsq, 0, 0.5, ylim=c(0,1))
plot(bexp, add=TRUE, col='red')
plot(blog, add=TRUE, col='blue')

c <- 0.7
b2sq <- Vectorize(bound2(lsquare, c))
b2exp <- Vectorize(bound2(lexp, c))
b2log <- Vectorize(bound2(llog, c))
bc <- B(lcost(c))		# The true regret of l_{0.5} (i.e., alpha)

est <- 0.4
plot(function(y) bc(y,est), 0, 1, ylim=c(0,1), lty=2)
plot(function(y) b2sq(y,est), add=TRUE)
plot(function(y) b2exp(y,est), add=TRUE, col='red')
plot(function(y) b2log(y,est), add=TRUE, col='blue')

plotBounds <- function(loss, c, est) {
	Bc <- B(lcost(c))
	bloss <- Vectorize(bound2(loss, c))
	plot(function(y) Bc(y,est), 0, 1, ylim=c(-1,1), lty=2)
	plot(function(y) psi(loss, c, Bc(y,est)), add=TRUE, col='red')
	plot(function(y) psi(loss, c, -Bc(y,est)), add=TRUE, col='red', lty=2)
}

