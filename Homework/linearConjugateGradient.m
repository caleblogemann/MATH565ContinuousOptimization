function [ x, k ] = linearConjugateGradient( A, b, x, TOL, MaxIter )
%LINEARCONJUGATEGRADIENT Summary of this function goes here
%   Detailed explanation goes here

    r = A*x - b;
    delta = r'*r;
    d = -r;
    k = 0;
    mstop = 1;
    while k < MaxIter && mstop
        k = k + 1;
        w = A*d;
        alpha = delta/(d'*w);
        x = x + alpha*d;
        r = r + alpha*w;
        deltaOld = delta;
        delta = r'*r;
        if sqrt(delta) < TOL
            mstop = 0;
        else
            beta = delta/deltaOld;
            d = -r + beta*d;
        end
    end
end

