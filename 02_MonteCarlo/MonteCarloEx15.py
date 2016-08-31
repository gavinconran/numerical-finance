###Exercise 15.
# Implement an Euler scheme for the Feller diffusion and plot the convergence
# of the price of the bond as the number of time steps become large.
# K = 2.1, theta = 0.09, v_0 = 0.09, sigma = 0.1, t = 2

### Does NOT include dW component in Euler_step ###

import numpy as np
import random
from matplotlib import pyplot

# function return analytic Bond Price
def analytic_B(K, sigma, theta, T, r_0):

    gamma = 0.5 * np.sqrt(K*K + 2*sigma*sigma)

    m_t = (2*K*theta/(sigma*sigma)) * np.log(gamma*np.exp(K*T/2.) / \
                                    (gamma*np.cosh(gamma*T) + K*np.sinh(gamma*T)/2.))

    n_t = np.sinh(gamma*T) / (gamma*np.cosh(gamma*T) + K*np.sinh(gamma*T)/2.) 
    
    return np.exp(m_t + n_t*r_0) 


random.seed(1234)

# function Euler step return the present solution plus a bit more
def euler_step(U, dt, theta, sigma, K):

    return U + dt*K*(theta-U) #+ sigma*np.sqrt(dt)*np.sqrt(abs(U)) *random.random()
    #return U + dt*K*(theta-U) + sigma*np.sqrt(dt)*np.sqrt(max(U, 0)) *random.random()
    #return U + dt*K*(theta-U) + sigma*np.sqrt(dt)*np.sqrt(U) *random.random()


# Construct Grid
T = 2.0
dt = 0.001
N = int(T/dt)+ 1

# Initialise solution array
t = np.linspace(0, T, N)
U = np.empty(N)
B_0 = 0.07
U[0] = B_0

# Set parameters
theta = 0.09
sigma = 0.1
K = 2.1

# time loop - Euler method
for n in range(0, N-1):
    U[n+1] = euler_step(U[n], dt, theta, sigma, K)

print("Approx. Price: ", U[-1])

# Compute Analytic Price
P_Analytic = analytic_B(K, sigma, theta,T, B_0)
print("Analytic Price: ", P_Analytic)

# Plot trajetory of Bond Prices over time
pyplot.figure(figsize=(8,6))
pyplot.grid(True)
pyplot.xlabel(r't', fontsize=18)
pyplot.ylabel(r'Bond Price', fontsize=18)
pyplot.title('Bond Price over time, delta t = %.3f' % dt, fontsize=18)
pyplot.plot(t,U, 'k-', lw=2);
pyplot.show()


# Convergence of Bond Price as time steps become large
N_List = range(102,5002, 10)
B = np.empty(len(N_List))
t = np.linspace(0, T, len(N_List))
index = 0          
for N in N_List:
    dt = T / (N-1)
    #print('dt: ', dt)
    U = np.empty(N)
    U[0] = 0.07

    # time loop - Euler method
    for n in range(0, N-1):
        U[n+1] = euler_step(U[n], dt, theta, sigma, K)

    B[index] = U[n+1]
    index += 1

# Plot convergence of Bond price as the number of time steps increases
pyplot.figure(figsize=(8,6))
pyplot.grid(True)
pyplot.xlabel(r'number of time steps', fontsize=18)
pyplot.ylabel(r'Bond Price', fontsize=18)
pyplot.title('Convergence of Bond Price', fontsize=18)
pyplot.plot(N_List,B, 'k-', lw=2);
pyplot.show()



    




