import numpy as np
def problem7(g, x0, TOL, MaxIter):
    x = np.zeros(MaxIter)
    x[0] = x0
    nIter = 0
    mstop = 1
    gx = g(x[nIter])
    delta = -np.square(gx)/(g(x[nIter] + gx) + gx)

    while nIter <= MaxIter and mstop:
        nIter+=1
        x[nIter] = x[nIter - 1] + delta
        gx = g(x[nIter])
        if abs(gx) < TOL and abs(delta) < TOL:
            mstop = 0
        else:
            delta = -np.square(gx)/(g(x[nIter] + gx) + gx)
    x = x[:nIter+1]
    return x
