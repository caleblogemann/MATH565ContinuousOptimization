def backtrace(f, gradf, x, p):
    a = 1
    rho = .8
    c = .1
    while f(x + a*p) > f(x) + c*a*np.dot(gradf(x), p):
        a = rho*a
    return a
