import numpy as np
def linearConjugateGradient(A, b, x, MaxIter, TOL):
    r = np.dot(A, x) - b
    delta = np.dot(r, r)
    p = -r
    k = 0
    mstop = 1
    while k < MaxIter and mstop:
        k+=1
        w = np.dot(A, p)
        alpha = delta/np.dot(p, w)
        x = x + alpha*p
        r = r + alpha*w
        deltaOld = delta
        delta = np.dot(r, r)
        if np.sqrt(delta) < TOL:
            mstop = 0
        else:
            beta = delta/deltaOld
            p = -r + beta*p
    return (x, k)
