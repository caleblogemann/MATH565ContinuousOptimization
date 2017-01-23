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

sol = problem7(g, h0, TOL, MaxIter)
print 'The depth of the sphere is {0:.10f}'.format(sol[-1])

x = np.linspace(1,10,10000)
plt.plot(x, g(x), '-', sol[-1], g(sol[-1]), 'o')

