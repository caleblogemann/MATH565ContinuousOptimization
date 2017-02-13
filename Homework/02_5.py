import numpy as np
import matplotlib.pyplot as plt
execfile('newton1d.py')
def steepestDescent(f, gradf, phid, phidd, x0, TOL, MaxIter):
    x = np.zeros([MaxIter+1, x0.size])
    x[0] = x0
    nIter = 0
    mstop = 1

    while nIter < MaxIter and mstop:
        p = -gradf(x[nIter])
        g = lambda a: phid(a, x[nIter], p)
        gd = lambda a: phidd(a, x[nIter], p)
        sol = newton1D(g, gd, 0, TOL, MaxIter)
        #if sol.size >= MaxIter:
            # print(nIter)
            # print(p)
            # print(x[nIter])
            # print(sol[-1])
            # raise Exception('Search for alpha did not converge')
        # alpha = backtrace(f, gradf, x[nIter], p)
        alpha = sol[-1]
        # print(nIter)
        # print(alpha)
        delta = alpha * p
        x[nIter+1] = x[nIter] + delta
        nIter+=1
        if sum(abs(delta)) < TOL and sum(abs(p)) < TOL:
            mstop = 0

    return x[:nIter+1]

def plotResults(f, sol, title):
    minX0 = .5
    maxX0 = 1.5
    minX1 = .5
    maxX1 = 1.5
    meshSize = 100

    x0list = np.linspace(minX0, maxX0, meshSize)
    x1list = np.linspace(minX1, maxX1, meshSize)
    X0, X1 = np.meshgrid(x0list, x1list)
    Z = np.zeros([meshSize,meshSize])
    for i in range(meshSize):
        for j in range(meshSize):
            Z[i,j] = f([X0[i,j], X1[i,j]])

    print('plotting')
    plt.figure()
    levels = [0.5, 2, 4, 8, 16, 32]
    contour = plt.contour(X0, X1, Z, levels, colors='k')
    plt.plot(sol[:,0], sol[:,1], '-k')
    plt.plot(sol[:,0], sol[:,1], 'ko')
    plt.clabel(contour, colors='k', fmt='%2.1f', fontsize=12)
    plt.xlabel('x1')
    plt.ylabel('x2')
    plt.title(title)
    plt.show()

f = lambda x: 20*np.square(x[1] - np.square(x[0])) + np.square(1 - x[0])
def gradf(x):
    dx0 = -80*x[0]*(x[1] - np.square(x[0])) - 2*(1 - x[0])
    dx1 = 40*(x[1] - np.square(x[0]))
    return np.array([dx0, dx1])

phid = lambda alpha, x, p: 40*(x[1] + alpha*p[1] - np.square(x[0] + alpha*p[0]))*(p[1] - 2*p[0]*(x[0] + alpha*p[0])) - 2*(1 - x[0] - alpha*p[0])*p[0]
phidd = lambda alpha, x, p: 40*np.square(p[1]) + 2*p[0]*(-80*p[1]*x[0] + p[0]*(1-40*x[1]-120*p[1]*alpha + 120*np.square(x[0] + p[0]*alpha)))

x0 = np.array([1.2, 1.2])
TOL = 1e-10
MaxIter = 1000
solSteepestDescent = steepestDescent(f, gradf, phid, phidd, x0, TOL, MaxIter)
plotResults(f, solSteepestDescent, 'Steepest Descent')
