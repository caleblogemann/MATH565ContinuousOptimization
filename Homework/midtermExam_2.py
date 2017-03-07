import numpy as np
execfile("trustRegionMethod.py")
execfile("doglegMethod.py")
f = lambda x: (x[0]**2 + x[1] - 11)**2 + (x[0] + x[1]**2 - 7)**2

def gradf(x):
    dx0 = 4*(x[0]**3 + x[0]*(x[1] - 11)) + 2*(x[0] + x[1]**2 - 7)
    dx1 = 2*(x[0]**2 + x[1] - 11) + 4*(x[1]**3 + x[1]*(x[0] - 7))
    return np.array([dx0, dx1])

def hessianf(x):
    dx0x0 = 12*x[0]**2 + 4*x[1] - 42
    dx0x1 = 4*x[0] + 4*x[1]
    dx1x1 = 12*x[1]**2 + 3*x[0] - 26
    return np.array([[dx0x0, dx0x1],[dx0x1, dx1x1]])

def plotResults(f, x0Range, x1Range, solList, title):
    minX0 = x0Range[0]
    maxX0 = x0Range[1]
    minX1 = x1Range[0]
    maxX1 = x1Range[1]
    meshSize = 100

    x0list = np.linspace(minX0, maxX0, meshSize)
    x1list = np.linspace(minX1, maxX1, meshSize)
    X0, X1 = np.meshgrid(x0list, x1list)
    Z = np.zeros([meshSize,meshSize])
    for i in range(meshSize):
        for j in range(meshSize):
            Z[i,j] = f([X0[i,j], X1[i,j]])

    plt.figure()
    levels = [10, 40, 80, 160, 320, 640]
    contour = plt.contour(X0, X1, Z, levels, colors='k')
    for i in range(len(solList)):
        plt.plot(solList[i][:,0], solList[i][:,1], '-k')
        plt.plot(solList[i][:,0], solList[i][:,1], 'ko')
    plt.clabel(contour, colors='k', fmt='%2.1f', fontsize=12)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(title)
    plt.show()

maxDelta = 2
Delta0 = 1

eta = .2
TOL = 1e-10
MaxIter = 100

x0 = [5.5, 5.5]
sol0 = trustRegionMethod(f, gradf, hessianf, doglegMethod, x0, maxDelta, Delta0, eta, TOL, MaxIter)
x0 = [-5.5, 5.5]
sol1 = trustRegionMethod(f, gradf, hessianf, doglegMethod, x0, maxDelta, Delta0, eta, TOL, MaxIter)
x0 = [5.5, -5.5]
sol2 = trustRegionMethod(f, gradf, hessianf, doglegMethod, x0, maxDelta, Delta0, eta, TOL, MaxIter)
x0 = [-5.5, -5.5]
sol3 = trustRegionMethod(f, gradf, hessianf, doglegMethod, x0, maxDelta, Delta0, eta, TOL, MaxIter)
solList = [sol0, sol1, sol2, sol3]

plotResults(f, [-6, 6], [-6, 6], solList, "")
