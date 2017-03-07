import numpy as np
import ipdb
def trustRegionMethod(f, gradf, B, minimizingFunction, x0, maxDelta, Delta0, eta, TOL, MaxIter):
    N = len(x0)
    x = np.zeros([MaxIter+1, N])
    x[0] = x0
    nIter = 0
    mstop = 1
    Delta = Delta0

    gradfx = gradf(x[nIter])
    while nIter < MaxIter and mstop:
        Bx = B(x[nIter])
        p = minimizingFunction(gradfx, Bx, Delta)

        m = lambda p: f(x[nIter]) + np.dot(p, gradfx) + (1.0/2.0)*np.dot(p, np.dot(Bx, p))
        deltaM = m(np.zeros(N)) - m(p)
        rho = (f(x[nIter]) - f(x[nIter] + p))/(deltaM)

        # modify trust region size
        if rho < 1.0/4.0:
            Delta = (1.0/4.0)*Delta
        elif rho > 3.0/4.0 and abs(np.linalg.norm(p) - Delta) < 1e-5:
            Delta = min(2*Delta, maxDelta)

        # Decide whether or not to reject step
        if rho > eta:
            x[nIter + 1] = x[nIter] + p
        else:
            x[nIter + 1] = x[nIter]

        nIter+=1
        if nIter % 100 == 0:
            print((nIter + 0.0)/MaxIter)
        gradfx = gradf(x[nIter])
        if np.linalg.norm(p) < TOL and np.linalg.norm(gradfx) < TOL:
            mstop = 0

    return x[:nIter+1]
