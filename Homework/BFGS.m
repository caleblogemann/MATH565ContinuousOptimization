function [] = BFGS(x0, H0, TOL, MaxIter)
    k = 0;
    mstop = 1;

    while k <= MaxIter && mstop
        p = -H;

        if
            mstop = 0;
        else
            k = k+1;
        end
    end

end
