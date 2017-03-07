execfile("trustRegionMethod.py")
execfile("doglegMethod.py")

d = np.array([[   0,  587, 1212,  701, 1936,  604,  748, 2139, 2182,  543,  762],
              [ 587,    0,  920,  940, 1745, 1188,  713, 1858, 1737,  597,  309],
              [1212,  920,    0,  879,  831, 1726, 1631,  949, 1021, 1494,  614],
              [ 701,  940,  879,    0, 1374,  968, 1420, 1645, 1891, 1220,  854],
              [1936, 1745,  831, 1374,    0, 2339, 2451,  347,  959, 2300, 1443],
              [ 604, 1188, 1726,  968, 2339,    0, 1092, 2594, 2734,  923, 1361],
              [ 748,  713, 1631, 1420, 2451, 1092,    0, 2571, 2408,  205, 1021],
              [2139, 1858,  949, 1645,  347, 2594, 2571,    0,  678, 2442, 1548],
              [2182, 1737, 1021, 1891,  959, 2734, 2408,  678,    0, 2329, 1451],
              [ 543,  597, 1494, 1220, 2300,  923,  205, 2442, 2329,    0,  898],
              [ 762,  309,  614,  854, 1443, 1361, 1021, 1548, 1451,  898,    0]])
N = len(d)
def f(v):
    ipdb.set_trace()
    x = v[::2]
    y = v[1::2]
    return sum([sum([((x[i] - x[j])**2 + (y[i] - y[j])**2 - d[i,j]**2)**2 for j in range(N)]) for i in range(N)])

def gradf(x):
    gradfx = np.zeros(2*N)
    for k in range(N):
        # \pd{f}{x_k}
        gradfx[2*k] = 8*sum([((x[2*k] - x[2*j])**2 + (x[2*k+1] - x[2*j+1])**2 - d[k,j]**2)*(x[2*k] - x[2*j]) for j in range(N)])
        # \pd{f}{y_k}
        gradfx[2*k+1] = 8*sum([((x[2*k] - x[2*j])**2 + (x[2*k+1] - x[2*j+1])**2 - d[k,j]**2)*(x[2*k+1] - x[2*j+1]) for j in range(N)])
    return gradfx

def hessianf(x):
    hessianfx = np.zeros([2*N, 2*N])
    for k in range(N):
        for i in range(N):
            if i == k:
                # \pd[2]{f}{x_k x_k}
                hessianfx[2*k, 2*k] = 8*sum([3*(x[2*k] - x[2*j])**2 + (x[2*k+1] - x[2*j+1])**2 - d[k,j]**2 for j in range(N)])
                # \pd[2]{f}{x_k y_k}
                hessianfx[2*k, 2*k+1] = 16*sum([(x[2*k] - x[2*j])*(x[2*k+1] - x[2*j+1]) for j in range(N)])
            else:
                # \pd[2]{f}{x_k x_i}
                hessianfx[2*k, 2*i] = -8*(3*(x[2*k] - x[2*i])**2 + (x[2*k+1] - x[2*i+1])**2 - d[k,i]**2)
                # \pd[2]{f}{x_k y_i}
                hessianfx[2*k, 2*i+1] = -16*(x[2*k+1] - x[2*i+1])*(x[2*k] - x[2*i])

        for i in range(N):
            if i == k:
                # \pd[2]{f}{y_k x_k}
                hessianfx[2*k+1, 2*k] = 16*sum([(x[2*k] - x[2*j])*(x[2*k+1] - x[2*j+1]) for j in range(N)])
                # \pd[2]{f}{y_k y_k}
                hessianfx[2*k+1, 2*k+1] = 8*sum([(x[2*k] - x[2*j])**2 + 3*(x[2*k+1] - x[2*j+1])**2 - d[k,j]**2 for j in range(N)])
            else:
                # \pd[2]{f}{y_k x_i}
                hessianfx[2*k+1, 2*i] = -16*(x[2*k+1] - x[2*i + 1])*(x[2*k] - x[2*i])
                # \pd[2]{f}{y_k y_i}
                hessianfx[2*k+1, 2*i+1] = -8*((x[2*k] - x[2*i])**2 + 3*(x[2*k+1] - x[2*i+1])**2 - d[k,i]**2)
    return hessianfx

def hessianf2(v):
    hessianfx = np.zeros([2*N, 2*N])
    x = v[::2]
    y = v[1::2]
    for m in range(2*N):
        k = m/2
        for n in range(2*N):
            i = n/2
            # derivative with respect to x_k
            if m % 2 == 0:
                # derivative with respect to x_i
                if n % 2 == 0:
                    # \pd[2]{f}{x_k}
                    if k == i:
                        hessianfx[m,n] = 8*sum([(3*(x[k] - x[j])**2 + (y[k] - y[j])**2 - d[k,j]**2) for j in range(N)])
                    # \mpd[2]{f}{\partial x_k \partial x_i}
                    else:
                        hessianfx[m,n] = -8*(3*(x[k] - x[i])**2 + (y[k] - y[i])**2 - d[k,i]**2)
                # derivative with respect to y_i
                else:
                    # \mpd[2]{f}{\partial x_k \partial y_k}
                    if k == i:
                        hessianfx[m, n] = 16*sum([(x[k] - x[j])*(y[k] - y[j]) for j in range(N)])
                    # \mpd[2]{f}{\partial x_k \partial y_i}
                    else:
                        hessianfx[m, n] = -16*(y[k] - y[i])*(x[k] - x[i])
            # derivative with respect to y_k
            else: 
                # derivative with respect to x_i
                if n % 2 == 0:
                    # \mpd[2]{f}{\partial y_k \partial x_k}
                    if k == i:
                        hessianfx[m, n] = 16*sum([(x[k] - x[j])*(y[k] - y[j]) for j in range(N)])
                    # \mpd[2]{f}{\partial y_k \partial x_i}
                    else:
                        hessianfx[m, n] = -16*(x[k] - x[i])*(y[k] - y[i])
                # derivative with respect to y_i
                else:
                    # \pd[2]{f}{y_k}
                    if k == i:
                        hessianfx[m, n] = 8*sum([((x[k] - x[j])**2 + 3*(y[k] - y[j])**2 - d[k,j]**2) for j in range(N)])
                    # \mpd[2]{f}{\partial y_k \partial y_i}
                    else:
                        hessianfx[m, n] = -8*((x[k] - x[i])**2 + 3*(y[k] - y[i])**2 - d[k,i]**2)
    return hessianfx

def plotResults(sol, title):
    plt.figure()
    for i in range(N):
        plt.plot(sol[:,2*i], sol[:,2*i+1], '-k')
        plt.plot(sol[:,2*i], sol[:,2*i+1], 'ko')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(title)
    plt.show()

maxDelta = 100
Delta0 = 100

eta = .2
TOL = 1e-10
MaxIter = 2000

x0 = np.zeros(2*N)
xCoordinates = 3000*np.random.rand(N) - 1500
yCoordinates = 1000*np.random.rand(N) - 500
x0[::2] = xCoordinates
x0[1::2] = yCoordinates

#sol = trustRegionMethod(f, gradf, hessianf, doglegMethod, x0, maxDelta, Delta0, eta, TOL, MaxIter)

plotResults(sol, "")
