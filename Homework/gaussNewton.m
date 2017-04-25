function [x, k] = gaussNewton(x0, J, r, TOL, MaxIter)
    two_n = length(x0);
    n = two_n/2;

    x = zeros(two_n, MaxIter+1);
    x(:,1) = x0;

    k = 1;
    mstop = 1;

    rx = r(x(:,k));
    Jx = J(x(:,k));
    while k <= MaxIter && mstop
        p = -(Jx'*Jx)\(Jx'*rx);
        x(:,k+1) = x(:,k) + p;

        k = k+1;
        rx = r(x(:,k));
        if norm(p) < TOL
            mstop = 0;
        else
            Jx = J(x(:,k));
        end
    end
    x = x(:,1:k);
end
