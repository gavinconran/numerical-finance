## NON-LINEAR optimisation (Lecture 8 - Kjell)
## Implied Volatility
# R function to compute Black-Scholes call price
bsc <- function(S, T, t, K, r, s, q) {
  d1 <- (log(S/K)+(r-q+0.5*s^2)*(T-t))/(s*sqrt(T-t))
  d2 <- d1-s*sqrt(T-t)
  S*exp(-q*(T-t))*pnorm(d1)-K*exp(-r*(T-t))*pnorm(d2)
}
# Since R treats all variables as vectors
bsc(50, 0.5, 0.0, 45, 0.06, 0.2, 0.02)
bsc(50, 0.5, 0.0, 45, 0.06, c(0.15, 0.2, 0.25), 0.02)

# Suppose the option sold for $7, find σ
# Plot f(σ) over a range of values and see where it crosses the x axis
sigmas <- seq(0.05, 0.5, by = 0.01)
fsig <- bsc(50, 0.5, 0.0, 45, 0.06, sigmas, 0.02) - 7
plot(sigmas, fsig, type = "l")
# The implied volatility is σ implied = 0.25

# To compute σ implied , had to evaluate Black-Scholes formula 46 times
# Computed answer still not very precise
bsc(50, 0.5, 0.0, 45, 0.06, 0.25, 0.02) - 7

# GOAL: compute σ implied to within a pre-specified tolerance with minimum number of function evaluations
# Methods are called NON-LINEAR SOLVERS

## Method 1: Bisection Method
# Function of ONE variable
# No error checking
bisection <- function(f, a, b, tol = 0.001) {
  while(b-a > tol) {
    c <- (a+b)/2
    if(sign(f(c)) == sign(f(a)))
      a <- c
    else
      b <- c
  }
  (a+b)/2
}

# Write f (σ) as a function of one variable
fsig <- function(sigma)
  bsc(50, 0.5, 0.0, 45, 0.06, sigma, 0.02) - 7

# Use bisection to solve f (σ) = 0
bisection(fsig, 0.1, 0.3)

# Check computed solution
bsc(50, 0.5, 0.0, 45, 0.06, 0.2511719, 0.02)

## Method 2; newton's Method for n-dimensions non-linear problem
# Let g(x, y) = 1 - (x - 10^4 - (y - 1)^4)

# Gradient of g(x, y) = F(x, y)
F <- function(x,y)
  c(4*(x-1)^3, 4*(y-1)^3)

# Gradient of F(x, y)
DF <- function(x,y)
  + diag(c(12*(x-1)^2, 12*(y-1)^2))

# Starting point
x <- c(0,0)

# Do 25 Newton iterations
for (i in 1:25)
  x <- x - solve(DF(x[1], x[2]), F(x[1], x[2]))

x
