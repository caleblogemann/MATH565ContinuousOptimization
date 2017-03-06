import numpy as np
import norms
def trustRegionMethod(f, gradf, B, minimizingFunction, x0, maxDelta, Delta0, eta, TOL, MaxIter):
    x = np.zeros(MaxIter+1)
    x[0] = x0
    nIter = 0
    mstop = 1

    gradfx = gradf(x[nIter]) 
    while nIter < MaxIter and mstop:
        Bx = B(x[nIter])
        p = minimizingFunction(gradfx, Bx, Delta)
        rho = (f(x[nIter]) - f(x[nIter] + p))/(-np.dot(gradfx, p) - (1/2)*np.dot(p, np.dot(Bx, p)))

        # modify trust region size
        if rho < 1/4:
            Delta = (1/4)*Delta
        elif rho > 3/4 and abs(np.linalg.norm(p) - Delta) < 1e-5:
            Delta = min(2*Delta, maxDelta)

        # Decide whether or not to reject step
        if rho > eta:
            x[nIter + 1] = x[nIter] + p
        else:
            x[nIter + 1] = x[nIter]

        nIter += 1
        gradfx = gradf(x[nIter])
        if np.linalg.norm(p) < TOL and np.linalg.norm(gradfx) < TOL:
            mstop = 0

    return x[:nIter+1]
