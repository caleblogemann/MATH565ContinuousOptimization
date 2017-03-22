TOL = 1e-6
MaxIter = 500
for n in [10, 20, 40, 80, 160, 320]
    x0 = zeros(n*n, 1);
    b = np.ones(n*n, 1);
    A = CreateA(n);
    [x, k] = linearConjugateGradient(A, b, x0, TOL, MaxIter)
    disp(k)
    ek = TOL
    e0 = sqrt((x - x0)*A*(x - x0))
    e = (ek/(2*e0))^(1.0/k)
    c = ((e + 1)/(1 - e))^2
    disp(c)
end
