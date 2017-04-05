function [xk] = secantMethod(f, x0, x1, TOL, MaxIter)
    xkm1 = x1;
    xkm2 = x0;
    fxkm1 = f(xkm1);
    fxkm2 = f(xkm2);
    k = 0;
    mstop = 1;
    while k < MaxIter && mstop
        k = k + 1;

        delta = -fxkm1*(xkm1 - xkm2)/(fxkm1 - fxkm2);
        xk = xkm1 + delta;

        xkm2 = xkm1;
        fxkm2 = fxkm1;
        xkm1 = xk;
        fxkm1 = f(xkm1);
        if norm(delta) < TOL && fxkm1 < TOL
            mstop = 0;
        end
    end
end
