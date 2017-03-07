import numpy as np
def cauchyPointMethod(g, B, Delta):
    gBg = np.dot(g, np.dot(B, g))
    if gBg <= 0:
        tau = 1
    else:
        tau = min(1, np.linalg.norm(g)**3/(Delta*gBg))
    return (-tau*Delta/np.linalg.norm(g))*g
