execfile("trustRegionMethod.py")
execfile("doglegMethod.py")

d = np.array([[   0,  587, 1212,  701, 1936,  604,  748, 2139, 2182,  543,  762],
              [ 587,    0,  920,  940, 1745, 1188,  713, 1858, 1737,  597,  309],
              [1212,  920,    0,  879,  831, 1726, 1631,  949, 1021, 1494,  614],
              [ 701,  940,  879,    0, 1374,  968, 1420, 1645, 1891, 1220,  854],
              [1936, 1745,  831, 1374,    0, 2339, 2451,  347,  959, 2300, 1443],
              [ 604, 1188, 1726,  968, 2339,    0, 1092, 2594, 2734,  923, 1361],
              [ 748,  713, 1631, 1420, 2451, 1092,    0, 2571, 2408,  205, 1021],
              [2139, 1858,  949, 1645,  357, 2594, 2571,    0,  678, 2442, 1548],
              [2182, 1737, 1021, 1891,  959, 2734, 2408,  678,    0, 2329, 1451],
              [ 543,  597, 1494, 1220, 2300,  923,  205, 2442, 2329,    0,  898],
              [ 762,  309,  614,  854, 1443, 1361, 1021, 1548, 1451,  898,    0]])
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
