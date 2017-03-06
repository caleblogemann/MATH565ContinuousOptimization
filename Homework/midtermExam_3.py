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

execfile("trustRegionMethod.py")
execfile("doglegMethod.py")

d = np.random.rand(5,5)
N = len(d)
def f(x):
    return sum([sum([((x[2*i] - x[2*j])**2 + (x[2*i+1] - x[2*j+1])**2 -d[i,j]**2)**2 for j in range(N)]) for i in range(N)])

def gradf(x):
    gradfx = np.zeros(2*N)
    for k in range(N):
        # \pd{f}{x_k}
        gradfx[2*k] = 8*sum([((x[2*k] - x[2*j])**2 + (x[2*k+1] - x[2*j+1])**2 - d[k,j])*(x[2*k] - x[2*j]) for j in range(N)])
        # \pd{f}{y_k}
        gradfx[2*k+1] = 8*sum([((x[2*k] - x[2*j])**2 + (x[2*k+1] - x[2*j+1])**2 - d[k,j])*(x[2*k+1] - x[2*j+1]) for j in range(N)])
    return gradfx

def hessianf(x):
    hessianfx = np.zeros([2*N, 2*N])
    for k in range(N):
        for i in range(N):
            if i == k:
                # \pd[2]{f}{x_k x_k}
                hessianfx[2*k, 2*k] = 8*sum([3*(x[2*k] - x[2*j])**2 + (x[2*k+1] - x[2*j+1])**2 - d[k,j] for j in range(N)])
                # \pd[2]{f}{x_k y_k}
                hessianfx[2*k, 2*k+1] = 16*sum([(x[2*k] - x[2*j])*(x[2*k+1] - x[2*j+1]) for j in range(N)])
            else:
                # \pd[2]{f}{x_k x_i}
                hessianfx[2*k, 2*i] = -8*(3*(x[2*k] - x[2*i])**2 + (x[2*k+1] - x[2*i+1])**2 - d[k,i])
                # \pd[2]{f}{x_k y_i}
                hessianfx[2*k, 2*i+1] = -16*(x[2*k+1] - x[2*i + 1])*(x[2*k] - x[2*i])

        for i in range(N):
            if i == k:
                # \pd[2]{f}{y_k x_k}
                hessianfx[2*k+1, 2*k] = 16*sum([(x[2*k] - x[2*j])*(x[2*k+1] - x[2*j+1]) for j in range(N)])
                # \pd[2]{f}{y_k y_k}
                hessianfx[2*k+1, 2*k+1] = 8*sum([(x[2*k] - x[2*j])**2 + 3*(x[2*k+1] - x[2*j+1])**2 - d[k,j] for j in range(N)])
            else:
                # \pd[2]{f}{y_k x_i}
                hessianfx[2*k+1, 2*i] = -16*(x[2*k+1] - x[2*i + 1])*(x[2*k] - x[2*i])
                # \pd[2]{f}{y_k y_i}
                hessianfx[2*k+1, 2*i+1] = -8*((x[2*k] - x[2*i])**2 + 3*(x[2*k+1] - x[2*i+1])**2 - d[k,i])
    return hessianfx

def plotResults(x0Range, x1Range, sol, title):
    plt.figure()
    for i in range(N):
        plt.plot(sol[:,2*i], sol[:,2*i+1], '-k')
        plt.plot(sol[:,2*i], sol[:,2*i+1], 'ko')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(title)
    plt.show()

maxDelta = 2
Delta0 = 1

eta = 1/8
TOL = 1e-10
MaxIter = 100

x0 = np.zeros(2*N)
xCoordinates = 3000*np.random.rand(N) - 1500
yCoordinates = 1000*np.random.rand(N) - 500
x0[::2] = xCoordinates
x0[1::2] = yCoordinates

plotResults(np.vstack((x0,x0)), "")
