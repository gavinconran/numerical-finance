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
    return U + dt*K*(theta-U) + sigma*np.sqrt(dt)*np.sqrt(abs(U)) *random.random() 


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

# Convergence of Bond Price as number of simulations become large
# Plot trajetory of Bond Prices over time (* num_sims)
pyplot.figure(figsize=(8,6))


num_sims = 50
sims = range(num_sims)
B_n = np.zeros(num_sims)
B_Price = B_0
for i in sims:
    # time loop - Euler method
    for n in range(0, N-1):
        U[n+1] = euler_step(U[n], dt, theta, sigma, K)
    B_Price += U[-1]
    pyplot.plot(t,U, 'k-', lw=2);
    B_n[i] = B_Price / (i+1)

pyplot.grid(True)
pyplot.xlabel(r't', fontsize=18)
pyplot.ylabel(r'Bond Price', fontsize=18)
pyplot.title('Bond Price over time, dt: %.3f, #sims: %d' % (dt, num_sims), fontsize=18)
pyplot.plot(t, U, 'k-', lw=2);
pyplot.show()


print("Approx. Price: ", B_n[-1])

# Compute Analytic Price
P_Analytic = analytic_B(K, sigma, theta,T, B_0)
print("Analytic Price: ", P_Analytic)


# Plot Convergence of Bond Prices wrt no. simulations
pyplot.figure(figsize=(8,6))
pyplot.grid(True)
pyplot.xlabel(r'# simulations', fontsize=18)
pyplot.ylabel(r'Bond Price', fontsize=18)
pyplot.title('Convergence of Bond Price wrt # simulations', fontsize=18)
sims = [x + 1 for x in sims]
pyplot.plot(sims, B_n, 'k-', lw=2);
pyplot.show()



# Convergence of Bond Price as number of time steps become large
N_List = range(2,5002, 1000)
B = np.empty(len(N_List))
t = np.linspace(0, T, len(N_List))

num_sims = 5
sims = range(num_sims)

sim = 0  
# for each # steps        
for N in N_List:
    dt = T / (N-1) # compute dt
    U = np.empty(N) # Initialise U (Euler array)
    U[0] = 0.07

    B_Price = 0
    for i in sims: # rum simulations
        # time loop - Euler method
        for n in range(0, N-1):
            U[n+1] = euler_step(U[n], dt, theta, sigma, K)
        B_Price += U[-1]
    B[sim] = B_Price / (i+1) # Expectation of Bond Price for this simulation
    sim += 1

# Plot convergence of Bond price as the number of time steps increases
pyplot.figure(figsize=(8,6))
pyplot.grid(True)
pyplot.xlabel(r'number of time steps', fontsize=18)
pyplot.ylabel(r'Bond Price', fontsize=18)
pyplot.title('Convergence of Bond Price wrt # time steps', fontsize=18)
pyplot.plot(N_List,B, 'k-', lw=2);
pyplot.show()




    



