## Langrange' Method + Newton's Method
# Example 1: Revisit the example from the Lagrange’s method slides
# max/min: 4x2 − 2x3
#  subject to: 2x1 − x2 − x3 = 0 AND x1^2 + x2^2 − 13 = 0

# R function for evaluating the function F
F <- function(x)
  c(-4+2*x[4]*x[1],
    6+2*x[4]*x[2],
    2*x[1]-x[2]-x[3],
    x[1]^2+x[2]^2-13)

# R function for evaluating the gradient of F
DF <- function(x)
  matrix(c(2*x[4], 0, 0, 2*x[1],
           0, 2*x[4], 0, 2*x[2], 
           2,-1, -1, 0, 
           2*x[1], 2*x[2], 0,  0),
         4, 4, byrow = TRUE)

# Starting point
x <- rep(1, 4)

#15 Newton iterations
for(i in 1:15)
  x <- x - solve(DF(x), F(x))
x

# Starting point (-1,-1,-1,-1)
x <- - rep(1, 4)

#15 Newton iterations
for(i in 1:15)
  x <- x - solve(DF(x), F(x))
x

# Lagrange + Newton's Method: Example 2
# minimize: 3x1 − 4x2 + x3 − 2x4
# subject to: −x2^2 + x3^2 + x4^2 = 1 AND 3x1^2 + x3^2 + 2x4^2 = 6

# Function to compute G(x , λ)
G <- function(x)
  c(3 + 6*x[6]*x[1], -4 - 2*x[5]*x[2],
    1 + 2*x[5]*x[3] + 2*x[6]*x[3],
    -2 + 2*x[5]*x[4] + 4*x[6]*x[4],
    -x[2]^2 + x[3]^2 + x[4]^2 - 1,
    3*x[1]^2 + x[3]^2 + 2*x[4]^2 - 6)

# Function to compute D G(x , λ)
DG <- function(x) {
  grad <- matrix(0, 6, 6)
  grad[1,] <- c(6*x[6], 0, 0, 0, 0, 6*x[1])
  grad[2,] <- c(0, -2*x[5], 0, 0, -2*x[2], 0)
  grad[3,] <- c(0, 0, 2*x[5] + 2*x[6], 0, 2*x[3], 2*x[3])
  grad[4,] <- c(0, 0, 0, 2*x[5] + 4*x[6], 2*x[4], 2*x[4])
  grad[5,] <- c(0, -2*x[2], 2*x[3], 2*x[4], 0, 0)
  grad[6,] <- c(6*x[1], 0, 2*x[3], 4*x[4], 0, 0)
  grad
}

# Starting point
x <- c(1, -1, 1, -1, 1, -1)
# Newton iterations
for(i in 1:25)
  x <- x - solve(DG(x), G(x))
# Numeric solution
x
# Does the point (x c , λ c ) correspond to a minimum or a maximum?
# Already know x c is a critical point of F (x , λ c )
# xc ∼ minimum if D2F (xc , λc ) positive definite (all matrix eigennvales > 0)
# xc ∼ maximum if D2F (xc , λc ) negative definite (all matrix eigennvales < 0)
# Already have D2F (xc , λc ), upper-left 4 × 4 block of DG(xc , λc )
round(DG(x), digits=3)
# The diagonal entries of the 4*4 matrix are the eigenvalues
# Follows that x c corresponds to a maximum

# Example 3: Maximum Expected Returns Optimization
# Vector of expected returns
mu = c(0.08, 0.10, 0.13, 0.15, 0.20)
# Asset returns covariance matrix
Sigma = matrix(c(0.019600, -0.007560, 0.012880, 0.008750, -0.009800, 
                 -0.007560, 0.032400, -0.004140, -0.009000, 0.0094506,
                 0.012880, -0.004140, 0.052900, 0.020125, 0.020125,
                 0.008750, -0.009000,  0.020125, 0.062500, -0.013125,
                 -0.009800, 0.009450, 0.020125, -0.013125, 0.122500),
               nrow=5, ncol=5)
# Target risk: σ P 2 = 0.25 2 = 0.0625
sigma2P = 0.25**2

# Function to compute G(w , λ)
G <- function(x, mu, Sigma, sigmaP2)
{
  n <- length(mu)
  c(mu + rep(x[n+1], n) + 2*x[n+2]*(Sigma %*% x[1:n]),
    sum(x[1:n]) - 1,
    t(x[1:n]) %*% Sigma %*% x[1:n] - sigmaP2)
}

# Function to compute D G(w , λ)
DG <- function(x, mu, Sigma, sigmaP2)
{
  n <- length(mu)
  grad <- matrix(0.0, n+2, n + 2)
  grad[1:n, 1:n] <- 2*x[n+2]*Sigma
  grad[1:n, n+1] <- 1
  grad[1:n, n+2] <- 2*(Sigma %*% x[1:n])
  grad[n+1, 1:n] <- 1
  grad[n+2, 1:n] <- 2*t(x[1:n]) %*% Sigma
  grad
}

# From starting point
x <- c(rep(0.5, 5), 1, 1)
x

# Newton iterations
for(i in 1:25)
  x <- x - solve(DG(x, mu, Sigma, 0.25^2),
                 G(x, mu, Sigma, 0.25^2))
x

# Recall: upper-left n × n block of D G(w , λ c ) ∼ Hessian of F (w , λ c )
DG(x, mu, Sigma, sigmaP2)[1:5, 1:5]

# Can check second order condition by computing eigenvalues
eigen(DG(x, mu, Sigma, sigmaP2)[1:5, 1:5])$values
# Eigenvalues are negative: Follows that x c corresponds to a maximum
# Computed w is a constrained maximum
t(x[1:5]) %*% mu

