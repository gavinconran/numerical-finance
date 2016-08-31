### https://en.wikipedia.org/wiki/Euler%E2%80%93Maruyama_method

import numpy as np
import matplotlib.pyplot as plt

num_sims = 5

y_init = 0.07
t_init = 0.0
t_end  = 2.

theta = 0.09
sigma = 0.1
K = 2.1

dt = 0.01
N = int((t_end-t_init)/dt)+ 1

dt   = float(t_end - t_init) / N
dW   = lambda dt: np.random.normal(loc = 0.0, scale = np.sqrt(dt))

t    = np.arange(t_init, t_end, dt)
y    = np.zeros(N)
y[0] = y_init

def euler_step(U, dt, theta, sigma, K):
    # abs(U) can be replace by max(U, 0): see page 52 of Antonie's notes
    return U + dt*K*(theta-U) + sigma*np.sqrt(dt)*np.sqrt(abs(U)) * np.random.normal(0, 1) 

for i_sim in range(num_sims):
    for i in xrange(1, t.size):
        y[i] = euler_step(y[i-1], dt, theta, sigma, K)
    plt.plot(t, y)

plt.show()
