import numpy as np
import matplotlib.pyplot as plt
execfile('01_7.py')

P = 2
a = 12.87
b = 0.1142
n = 1
R = 0.08205
T = 313
g = lambda V: (P + a/np.square(V))*(V - b) - n*R*T
V0 = 12

TOL = 1e-10
MaxIter = 1000

sol = problem7(g, V0, TOL, MaxIter)
print 'The volume of the isobutane is {0:.10f}'.format(sol[-1])

x = np.linspace(1,20,10000)
plt.plot(x, g(x), '-', sol[-1], g(sol[-1]), 'o')

