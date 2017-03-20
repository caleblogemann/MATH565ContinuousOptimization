import numpy as np
from scipy.linalg import hilbert
execfile('linearConjugateGradient.py')

TOL = 1e-6
MaxIter = 100
for n in [5, 8, 12, 20]:
    A = hilbert(n)
    b = np.ones(n)
    x = np.zeros(n)
    (x, k) = linearConjugateGradient(A, b, x, MaxIter, TOL)
    print(k)

