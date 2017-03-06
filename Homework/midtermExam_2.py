import trustRegionMethod
import doglegMethod

f = lambda x: (x[0]^2 + x[1] - 11)**2 + (x[0] + x[1]**2 - 7)**2

def gradf(x):
    dx0 = 4*(x[0]**3 + x[0]*(x[1] - 11)) + 2*(x[0] + x[1]**2 - 7)
    dx1 = 2*(x[0]**2 + x[1] - 11) + 4*(x[1]*(x[0]-7) + x[1]**3)
    return np.array([dfdx0, dfdx1])

def hessianf(x):
    dx0x0 = 12*x[0]**2 + 4*x[1] - 42
    dx0x1 = 4*x[0] + 4*x[1]
    dx1x1 = 12*x[1]**2 + 4*x[0] - 26
    return np.array([[dx0x0, dx0x1], [dx0x1, dx1x1]])

maxDelta = 
Delta0 = 
eta = 1/8
TOL = 1e-10
MaxIter = 1e3

x0 = 

sol = trustRegionMethod(f, gradf, hessianf, doglegMethod, x0, maxDelta, Delta0, eta, TOL, MaxIter)

