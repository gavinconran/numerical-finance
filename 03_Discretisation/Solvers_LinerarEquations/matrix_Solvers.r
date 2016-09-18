## Examples from matrixReviewPowerPoint.pdf
#matrix
matA = matrix(data=c(1,2,3,4,5,6), nrow=2, ncol = 3)
class(matA)
dim(matA)
nrow(matA)

ncol(matB)

# Vector
xvec = c(1,2,3)
class(xvec)

## Coercion vector --> matrix
xvec.mat = as.matrix(xvec)
xvec.mat

class(xvec.mat)

# Check for Symmetry
matS = matrix(c(1,2,2,1),2,2)
matS == t(matS)

# Matrix aarithemtic operators (+, -, %*% (dot product), scalar multipication, * (entry wise multiplication))
matA = matrix(c(4,9,2,1),2,2,byrow=T)
matB = matrix(c(2,0,0,7),2,2,byrow=T)
# Addition
matC = matA + matB

# Subtraction
matC = matA - matB

# Multiplication %*% --> dot product
matA = matrix(1:4, 2,2, byrow=T)
matB = matrix(5:8, 2,2,byrow=T)
matC = matA%*%matB
matC

# Scalar multiplication
2 * matA

# Entry wise multiplication
matA * matB

# IDENTITY MATRIX
# create identity matrix
matI = diag(2)
matI

matI %*% matA
matA %*%matI

# dot  product of two vectors
vecA = c(1,8)
vecB = c(2, -7)
vecA%*%vecB

# MATRIX INVERSION & Transpose
# Transpose
t(matA)
t(xvec.mat)

# Inversion
# inversion is a little more complex, partly because the function you’d want to use has a non-obvious name: solve.
# The reason that solve is called solve is that it’s a general purpose function you can use to solve matrix equations without 
# wasting time computing the full inverse, which is often inefficient. 
matA.inv = solve(matA)
matA.inv

# definition of a matrix’s inverse is that the product of the matrix and its inverse is the identity matrix, if the inverse exists
matA%*%matA.inv # returns the Identity matrix
matA.inv%*%matA # also returns the Identity matrix

## Solving systems of equations
A <- matrix(c(2, 4, -2, 4, 9, -3, -2, -3, 7), 3, 3)
b <- c(2, 8, 10)
solve(A,b)

## LU decomposition
library(matrixcalc)
A <- matrix(c(2,3,1,4,7,5,0,-2,2), ncol=3, nrow=3,byrow=T)
luA <- lu.decomposition( A )
L <- luA$L
U <- luA$U

### MATRIX FACTORISATION
## SINGULAR VALUE FACTORISATION: Singular value Decomposition (SVD)
# Eaxmple 1: data
hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, "+") }
library(MASS)
X <- hilbert(9)[, 1:6]

svdR <- svd(X)

U <- svdR$u
S <- diag(svdR$d)
V <- svdR$v

all.equal(U %*% S %*% t(V), X)

# Example 2: Image (tux)
library(imager)
im <- load.image('Tux.jpg')
layout(t(1:2)) 
plot(im,main="Original image")
im.bw <- grayscale(im)
plot(im.bw,main="BW Image image")

svdIM <- svd(im.bw)
d <- diag(svdIM$d)
dim(d)

u <- svdIM$u
v <- svdIM$v
layout(t(1:1)) 
plot(1:length(svdIM$d), svdIM$d)


## EIGENVALUES & EIGENVECTORS
# Eigenvectors are surely the bane of every starting student of linear algebra, 
# though their considerable power to simplify problems makes them the darling of every applied mathematician. 
# Thankfully, R makes it easy to get these for every matrix:
eig <- eigen(matA)
S <- eig$vectors
lambda <- eig$values

## QR FACTORISATION: Least Squares Example
# First, get the data
library(quantmod)
getSymbols(c("C", "GSPC"))
citi <- c(coredata(monthlyReturn(C["2010"])))
sp500 <- c(coredata(monthlyReturn(GSPC["2010"])))

# The x variable is sp500, bind a column of ones to get matrix X
X <- cbind(1, sp500)

# Compute QR factorization of X and extract the Q and R matrices
qrX <- qr(X)
Q <- qr.Q(qrX, complete = TRUE)
R <- qr.R(qrX, complete = TRUE)

# Compute u = Q T y
u <- t(Q) %*% citi

# Solve for αˆ and β ˆ
backsolve(R[1:2,1:2], u[1:2])

# Compare with built-in least squares fitting function
coef(lsfit(sp500, citi))
