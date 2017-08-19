
from pylab import *

# read in three columns from file and unpack into 3 arrays:
n,int_approx,error = loadtxt('mc_quad_error.txt',unpack=True)

figure(1)
clf()
loglog(n,error,'-o',label='Monte-Carlo')
loglog([1,1e7],[1,sqrt(1e-7)],'k',label='1 / sqrt(N)')
legend()
xlabel('number of MC points used')
ylabel('abs(error)')
title('Log-log plot of relative error in MC quadrature')
savefig('mc_quad_error.png')
