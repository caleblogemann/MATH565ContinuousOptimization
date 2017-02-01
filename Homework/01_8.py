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

TOL = 1e-10
MaxIter = 1000

V0 = 12
sol1 = problem7(g, V0, TOL, MaxIter)
V0 = .175
sol2 = problem7(g, V0, TOL, MaxIter)
V0 = .3375
sol3 = problem7(g, V0, TOL, MaxIter)

print 'The possible volumes of the isobutane are {0:.10f}, {1:.10f}, and {2:.10f}'.format(sol1[-1], sol2[-1], sol3[-1])
x = np.linspace(.15,15,10000)
plt.plot(x, g(x), '-', [sol1[-1], sol2[-1], sol3[-1]], [g(sol1[-1]), g(sol2[-1]), g(sol3[-1])], 'o')

