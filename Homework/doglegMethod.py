import numpy as np
def doglegMethod(g, B, Delta):
    # exact minimizer of approximation function, m
    pB = -np.linalg.solve(B, g)
    if np.linalg.norm(pB) < Delta:
        # print('Exact')
        return pB

    # steepest descent minimizer
    pU = - np.dot(g, g)/np.dot(g, np.dot(B, g)) * g
    if np.linalg.norm(pU) > Delta:
        # print('Steepest Descent')
        return Delta/np.linalg.norm(pU) * pU

    print('Dogleg Path')
    # use dogleg path
    # solve (pU + alpha(pB - pU))^T (pU + alpha(pB - pU)) = Delta^2
    a = np.dot(pB - pU, pB - pU)
    b = 2*np.dot(pU, pB - pU)
    c = np.dot(pU, pU) - Delta**2

    alpha0 = (-b + sqrt(b**2 - 4*a*c))/(2*a)
    alpha1 = (-b - sqrt(b**2 - 4*a*c))/(2*a)

    if alpha0 > 0 and alpha0 < 1:
        return pU + alpha1*(pB - pU)
    else:
        return pU + alpha0*(pB - pU)
