#########Running an explicit finite-difference scheme for the heat equation with boundary conditions

from math import *
from numpy import *
from pylab import *


def BCSpace(x): ## Space-boundary conditions
	if x<=0.5 and x>= 0.0:
		return x*x
	elif x>0.5 and x<=1.0:
		return (1.0 - x)*(1.0 - x)
	else:
		return 0.0
		
def BCTime(x): ## Time-boundary conditions
	return 0.0
	
def plotInit(X,U): ## Plot any intermediate solution
	p = size(U)
	Utemp = [BCSpace(0.0)]
	for i in range(p):
		Utemp.append(U[i])
	Utemp.append(BCSpace(xu))
	return plot(X, Utemp, 'b', label="initial boundary condition")

if __name__ == "__main__":
	sigma = sqrt(0.2)
	xl = 0.0
	xu = 1.0
	T = 1.0
	m = 100 ## discretisation in space
	p = m-1
	n = 10000 ## discretisation in time 
	dt = T/n
	dx = (xu - xl) / m
	alpha = dt / (dx*dx)
	print "dx = ", dx
	print "dt = ", dt
	
	A = zeros((p,p))
	A[0,0] = 1.0 - alpha*sigma*sigma ## diagonal term
	A[0, 1] = 0.5*alpha*sigma*sigma ## upper diagonal term
	
	for i in range(1,p-1):
		A[i, i] = 1.0 - alpha*sigma*sigma ## diagonal term
		A[i, i-1] = 0.5*alpha*sigma*sigma ## lower diagonal term
		A[i, i+1] = 0.5*alpha*sigma*sigma ## upper diagonal term
	
	A[p-1,p-1] = 1.0 - alpha*sigma*sigma ## diagonal term
	A[p-1, p-2] = 0.5*alpha*sigma*sigma ## lowerdiagonal term
	
	X = zeros(m+1)
	B = zeros(p)
	B[0] = BCSpace(xl)
	B[p-1] = BCSpace(xu)
	U = zeros(p)
	
	X[0] = xl
	X[m] = xu
	for i in range(p): ## initialise the vector u at time 0
		X[i+1] = xl + (i+1.0)*dx
		U[i] = BCSpace(X[i+1])
	
	## Plot initial boundary condition
	pinit = plotInit(X, U)
	pinit
	#show()


	Output = [BCSpace(0.0)]
	
	for i in range(n):
		U = add(A.dot(U), 0.5*alpha*sigma*sigma*B)
		#plotInit(X, U)
		#print i, " / ",  A
		
	for i in range(p):
		Output.append(U[i])
	Output.append(BCSpace(xu))
	pfinal = plot(X, Output, 'r', label="final solution")
	legend(loc=1)
	title("Explicit finite difference scheme for the heat equation");
	show()
