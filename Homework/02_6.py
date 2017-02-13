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

import numpy as np
import matplotlib.pyplot as plt
def newtonMultiD(f, gradf, hessianf, x0, TOL, MaxIter):
    x = np.zeros([MaxIter+1, x0.size])
    x[0] = x0
    nIter = 0
    mstop = 1
    gradfx = gradf(x0)
    while nIter < MaxIter and mstop:
        p = -np.dot(np.linalg.inv(hessianf(x[nIter])), gradfx)
        x[nIter+1] = x[nIter] + p
        nIter+=1
        gradfx = gradf(x[nIter])
        if sum(np.square(p)) < TOL and sum(np.square(gradfx)) < TOL:
            mstop = 0

    return x[:nIter+1]

f = lambda x: 20*np.square(x[1] - np.square(x[0])) + np.square(1 - x[0])
def gradf(x):
    dx0 = -80*x[0]*(x[1] - np.square(x[0])) - 2*(1 - x[0])
    dx1 = 40*(x[1] - np.square(x[0]))
    return np.array([dx0, dx1])

def hessianf(x):
    dx0x0 = 2 + 240*np.square(x[0]) - 80*x[1]
    dx0x1 = -80*x[0]
    dx1x1 = 40
    return np.array([[dx0x0, dx0x1], [dx0x1, dx1x1]])

x0 = np.array([1.2, 1.2])
TOL = 1e-10
MaxIter = 1000
solNewton = newtonMultiD(f, gradf, hessianf, x0, TOL, MaxIter)
plotResults(f, solNewton, 'Newton\'s Method')
