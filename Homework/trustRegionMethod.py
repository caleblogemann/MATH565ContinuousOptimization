import numpy as np
import norms
def trustRegionMethod(f, gradf, mk, x0, maxDelta, Delta0, eta, TOL, MaxIter):
    x = np.zeros(MaxIter+1)
    x[0] = x0
    nIter = 0
    mstop = 1

    while nIter < MaxIter and mstop:
        m = lambda 0: mk()
        p = 
        rho = (f(x[nIter]) - f(x[nIter] + p))/(m(0) - m(p))


        # modify trust region size
        if rho < 1/4:
            Delta = (1/4)*Delta
        elif rho > 3/4 and abs(norm2(p) - Delta) < 1e-5:
            Delta = min(2*Delta, maxDelta)

        # Decide whether or not to reject step
        if rho > eta:
            x[nIter + 1] = x[nIter] + p
        else:
            x[nIter + 1] = x[nIter]


        if sum(np.square(p)) < TOL and sum(np.square(gradfx)) < TOL:
            mstop = 0

    return x[:nIter+1]
