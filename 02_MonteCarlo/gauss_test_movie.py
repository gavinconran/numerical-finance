# Box Muller transform 
# Generates standard iid random numbers
# from uniformly distributed random numbers

import random, math, pylab

def gauss_test(sigma):
    # flat input distributions of phi and Upsilon
    # get transformed into x and y which are Gaussians and independent (iid)
    # This show that the substitution rules (the changes of variables in integrals)
    # also apply to the samples.
    phi = random.uniform(0.0, 2.0 * math.pi)
    Upsilon = random.uniform(0.0, 1.0)
    # intermediates variables
    Psi = - math.log(Upsilon)
    r = sigma * math.sqrt(2.0 * Psi)
    # and finally the gaussian sampled x and y numbers
    x = r * math.cos(phi)
    y = r * math.sin(phi)
    return [x, y]
    #return [r] # Distribution of r looks a bit like Boltzmann distribution

# exact distrubution:
list_x = [i * 0.1 for i in xrange(-40, 40)]
list_y = [math.exp(- x ** 2 / 2.0) / (math.sqrt(2.0 * math.pi)) for x in list_x]
# sampled distribution:
n_sampled_pairs = 50000
data = []
for sample in xrange(n_sampled_pairs):
        data += gauss_test(1.0)
# graphics output
pylab.plot(list_x, list_y, color='k', label='exact')
pylab.hist(data, bins=150, normed=True, color='r', histtype='step', label='sampled')
pylab.legend()
pylab.title('Sampling of the gaussian distribution\n(gauss_test_movie.py)')
pylab.xlabel('$x$', fontsize=14)
pylab.ylabel('$\pi(x)$', fontsize=14)
pylab.savefig('plot-gauss_test.png')

