function [x, k] = BFGS(x0, gradf, H0, alphaSearch, TOL, MaxIter)
    n = length(x0);
    k = 1;
    mstop = 1;

    H = H0;
    x = zeros(n, MaxIter+1);
    x(:,1) = x0;

    gradfx = gradf(x(:,k));

    I = eye(n);

    while k <= MaxIter && mstop
        p = -H*gradfx;
        alpha = alphaSearch(x(:,k), p);
        s = alpha*p;
        x(:,k+1) = x(:,k) + s;

        k = k+1;
        gradfxOld = gradfx;
        gradfx = gradf(x(:,k));
        if norm(gradfx) < TOL
            mstop = 0;
        else
            y = gradfx - gradfxOld;
            rho = 1/(y'*s);
            H = (I - rho*s*y')*H*(I - rho*y*s') + rho*s*s';
        end
    end
    x = x(:,1:k);
end
