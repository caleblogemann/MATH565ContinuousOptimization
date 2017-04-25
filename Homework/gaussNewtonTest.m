n = 7;
x = [0.038, 0.194, 0.425, 0.626, 1.253, 2.500, 3.740]';
y = [0.050,	0.127, 0.094, 0.2122, 0.2729, 0.2665, 0.3317]';
rj = @(j, b) y(j) - b(1)*x(j)/(b(2) + x(j));
r = @(b) arrayfun(@(j) rj(j,b),1:n)';
J = @(b) [arrayfun(@(j) -x(j)/(b(2)+x(j)), 1:n)', arrayfun(@(j) (b(1)*x(j))/(b(2)+x(j))^2, 1:n)'];
beta0 = [.9, .2];

sol = gaussNewton(beta0, J, r, TOL, MaxIter);
