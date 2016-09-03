###Exercise 15.
# Implement an Euler scheme for the Feller diffusion and plot the convergence
# of the price of the bond as the number of time steps become large.
# K = 2.1, theta = 0.09, v_0 = 0.09, sigma = 0.1, t = 2

### Include dW component in Euler_step ###
### got some help from https://en.wikipedia.org/wiki/Euler%E2%80%93Maruyama_method

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


# function Euler step return the present solution plus a bit more
def euler_step(U, dt, theta, sigma, K):
    # abs(U) can be replace by max(U, 0): see page 52 of Antonie's notes
    return U + dt*K*(theta-U) + sigma*np.sqrt(dt)*np.sqrt(abs(U)) * np.random.normal(0, 1) 

# Construct Grid
T = 2.0
dt = 0.01
N = int(T/dt)+ 1

# Initialise solution array
t = np.linspace(0, T, N)
U = np.empty(N) # solution array
r_0 = 0.07 # initial interest rate
U[0] = r_0 # Initial value 

# Set parameters
theta = 0.09 # long term mean of interest rate
sigma = 0.1  # variance of interest rate
K = 2.1      # mean-reversion strength

### Simulation of Stochastic Interest rates
# Plot trajetory of Stochastic interest rates over time (* num_sims)
pyplot.figure(figsize=(8,6))

num_sims = 5 # number of simulations
sims = range(num_sims)  # list of simulation runs
r_n = np.zeros(num_sims)  #  computed rate per simulation run
for i in sims:
    # time loop - Euler method
    for n in range(0, N-1):
        U[n+1] = euler_step(U[n], dt, theta, sigma, K)
    pyplot.plot(t,U, '-', lw=1);

pyplot.grid(True)
pyplot.xlabel(r't', fontsize=18)
pyplot.ylabel(r'Rate', fontsize=18)
pyplot.title('Rate over time, dt: %.3f, #sims: %d' % (dt, num_sims), fontsize=18)
pyplot.plot(t, U, 'k-', lw=2);

### Convergence of Bond Price as number of simulations become large
num_sims = 1000
sims = range(num_sims)  # range of number of runs per simulation
r_n = np.zeros(num_sims) # computed rate per simulation run
r_n_Totals = np.zeros(num_sims)
r_Total = 0
for i in sims:
    # time loop - Euler method
    for n in range(0, N-1):
        U[n+1] = euler_step(U[n], dt, theta, sigma, K)
    r_Total += U[-1]
    r_n_Totals[i] = r_Total / (i+1)
    r_n[i] = U[-1] # Record the rate for this simulation

r_Expectation = np.mean(r_n) # compute Expected value for rate
P_Expectation = np.exp(-r_Expectation) # Compute the Expected vale of the Bond Price from mean of the simulated rates


# Compute Analytic Price
P_Analytic = analytic_B(K, sigma, theta,T, r_0)
print("Analytic Price: ", P_Analytic)
print("Approx. Price: ", P_Expectation)

# Plot Histogram of Simulated Interest Rates at time T
pyplot.figure(figsize=(8,6))
pyplot.ylabel(r'Number of Hits', fontsize=18)
pyplot.xlabel(r'Interest Rates', fontsize=18)
pyplot.title('Histogram of Simulated Interest Rates. # sims: %d' % num_sims, fontsize=18)
pyplot.hist(r_n, 50)

# Plot Convergence of Bond Prices wrt no. simulations
pyplot.figure(figsize=(8,6))
pyplot.grid(True)
pyplot.xlabel(r'# simulations', fontsize=18)
pyplot.ylabel(r'Bond Price', fontsize=18)
pyplot.title('Convergence of Bond Price wrt # simulations', fontsize=18)
sims = [x + 1 for x in sims] # shift from 0 -> n-1 to  1 -> n (cosmetic)
B_n = [np.exp(-r) for r in r_n_Totals] # Compute Bond prices from simulated rates
pyplot.plot(sims, B_n, 'k-', lw=2, label='Approx. Price ');
pyplot.plot([1, num_sims], [P_Analytic, P_Analytic], 'r-', label='Analytic Price ')
pyplot.legend(loc='lower right');
pyplot.ylim(0.85, 0.95)


### Convergence of Bond Price as number of time steps become large
N_List = range(2, 22) #102)
r_n = np.empty(len(N_List)) 
t = np.linspace(0, T, len(N_List))

num_sims = 1000
sims = range(num_sims)

sim = 0  
# for each # steps        
for N in N_List:
    dt = T / (N-1) # compute dt
    U = np.empty(N) # Initialise U (Euler array)
    U[0] = r_0

    r_Total = 0 
    for i in sims: # rum simulations
        # time loop - Euler method
        for n in range(0, N-1):
            U[n+1] = euler_step(U[n], dt, theta, sigma, K)
        r_Total += U[-1]
    r_n[sim] = r_Total / (i+1) # Expectation of Bond Price for this simulation
    sim += 1

# Plot convergence of Bond price as the number of time steps increases
pyplot.figure(figsize=(8,6))
pyplot.grid(True)
pyplot.xlabel(r'number of time steps', fontsize=18)
pyplot.ylabel(r'Bond Price', fontsize=18)
pyplot.title('Convergence of Bond Price wrt # time steps', fontsize=18)
B_n = [np.exp(-r) for r in r_n] # Compute Bond prices from simulated rates
pyplot.plot(N_List, B_n, 'k-', lw=2, label='Approx. Price ');
pyplot.plot([0, max(N_List)], [P_Analytic, P_Analytic], 'r-', label='Analytic Price ')
pyplot.legend(loc='lower right');
pyplot.show()




    




