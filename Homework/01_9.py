import numpy as np
import matplotlib.pyplot as plt
execfile('01_7.py')

r = 5.0
rho = 1.0
sigma = .6
g = lambda h: 1.0/3.0 * math.pi * (3.0*r*np.square(h) - np.power(h, 3.0))*rho - 4.0/3.0*math.pi*np.power(r, 3.0)*sigma
h0 = 5.67

TOL = 1e-10
MaxIter = 1000

h0 = 5.67
sol1 = problem7(g, h0, TOL, MaxIter)
h0 = -4
sol2 = problem7(g, h0, TOL, MaxIter)
h0 = 13
sol3 = problem7(g, h0, TOL, MaxIter)

print 'The possible depths of the sphere are {0:.10f}, {1:.10f}, and {2:.10f}'.format(sol1[-1], sol2[-1], sol3[-1])

x = np.linspace(-5,15,10000)
plt.plot(x, g(x), '-', [sol1[-1], sol2[-1], sol3[-1]], [g(sol1[-1]), g(sol2[-1]), g(sol3[-1])], 'o')

