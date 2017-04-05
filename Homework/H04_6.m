n = 2;
f = @(x) (x(1)^2 + x(2) - 11)^2 + (x(1) + x(2)^2 - 7)^2;
gradf = @(x) [2*(x(1)^2 + x(2) - 11)*2*x(1) + 2*(x(1) + x(2)^2 - 7); 2*(x(1)^2 + x(2) - 11) + 2*(x(1) + x(2)^2 - 7)*2*x(2)];
H0 = eye(2);
TOL = 1e-5;
MaxIter = 1000;
phi = @(a, x, p) f(x + a*p);
phid = @(a, x, p) gradf(x + a*p)'*p;
alphaSearch = @(x, p) secantMethod(@(a) phid(a, x, p), 0, .1, TOL, MaxIter);

x01 = [1; 1];
[x1, k1] = BFGS(x01, gradf, H0, alphaSearch, TOL, MaxIter);
x02 = [-2; 5];
[x2, k2] = BFGS(x02, gradf, H0, alphaSearch, TOL, MaxIter);
x03 = [-2; -2];
[x3, k3] = BFGS(x03, gradf, H0, alphaSearch, TOL, MaxIter);
x04 = [5; -5];
[x4, k4] = BFGS(x04, gradf, H0, alphaSearch, TOL, MaxIter);

X = linspace(-6, 6, 1000);
Z = zeros(1000,1000);
for i=1:1000
    for j=1:1000
        Z(j,i) = f([X(i); X(j)]);
    end
end
hold on
contourf(X, X, Z, 30)
plot(x1(1,:), x1(2,:),'k', 'LineWidth', 5)
plot(x2(1,:), x2(2,:),'k', 'LineWidth', 5)
plot(x3(1,:), x3(2,:),'k', 'LineWidth', 5)
plot(x4(1,:), x4(2,:),'k', 'LineWidth', 5)
hold off
