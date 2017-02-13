import numpy as np
def newton1D(g, gd, x0, TOL, MaxIter):
    x = np.zeros(MaxIter+1)
    x[0] = x0
    nIter = 0
    mstop = 1

    gx = g(x[nIter])
    while nIter < MaxIter and mstop:
        gdx = gd(x[nIter])
        delta = -gx/gdx
        x[nIter+1] = x[nIter] + delta
        nIter+=1
        gx = g(x[nIter])
        if abs(gx) < TOL and abs(delta) < TOL:
            mstop = 0
    return x[:nIter + 1]
