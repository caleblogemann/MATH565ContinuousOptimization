from scipy.sparse import spdiags
from scipy.sparse import lil_matrix
from scipy.sparse import identity

execfile('linearConjugateGradient.py')
import ipdb

def CreateA(n1):
    e1= np.ones(n1)
    e2= np.zeros(n1)
    T = spdiags([e1, -4*e1, e1], [-1, 0, 1], n1, n1)
    I = identity(n1)

    A = lil_matrix((n1*n1,n1*n1))
    for i in range(n1):
        ma = i*n1
        mb = (i+1)*n1
        A[ma:mb,ma:mb] = T

    for i in range(n1-1):
        ma = i*n1
        mb = (i+1)*n1
        mc = (i+1)*n1
        md = (i+2)*n1
        A[ma:mb,mc:md] = I
        A[mc:md,ma:mb] = I

    A = -A
    return(A)

TOL = 1e-6
MaxIter = 500
for n in [10, 20, 40, 80]:
    x0 = np.zeros(n*n)
    b = np.ones(n*n)
    A = CreateA(n).toarray()
    (x, k) = linearConjugateGradient(A, b, x0, MaxIter, TOL)
    print(k)
    ek = TOL
    e0 = np.linalg.norm(x - x0)
    e = (ek/(2*e0))**(1.0/k)
    cond = ((e + 1)/(1 - e))**2
    print(cond)
