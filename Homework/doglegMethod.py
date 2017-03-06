import numpy as np
def doglegMethod(g, B, Delta):

    # exact minimum of m the approximation of f
    pB = -np.linalg.solve(B, g)
    # if exact solution is inside trust region use that
    if np.linalg.norm(pB) < Delta:
        return pB

    # steepest descent solution
    pU = -np.dot(g, g)/np.dot(g, np.dot(B, g)) * g
    # if steepest descent solution is outside trust region
    # use intersection of trust region boundary and steepest descent direction
    if np.linalg.norm(pU2) > Delta:
        return Delta*pU/norm(pU)

    # if steepest descent solution is inside trust region use dogleg path
    # solve ||pU + alpha(pB - pU)||^2 = Delta, which is quadractic in alpha
    # alpha needs to be between 0 and 1
    a = np.dot(pB - pU, pB - pU)
    b = 2*np.dot(pU, pB - pU)
    c = np.dot(pU, pU) - Delta

    alpha0 = -b + sqrt(b**2 - 4*a*c)/(2*a)
    alpha1 = -b - sqrt(b**2 - 4*a*c)/(2*a)

    if alpha0 > 0 and alpha0 < 1:
        return pU + alpha0*(pB - pU)
    else:
        return pU + alpha1(pB - pU)
