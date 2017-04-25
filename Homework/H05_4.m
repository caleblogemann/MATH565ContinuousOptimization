n = 2;
t = [0.0, 0.3, 0.6, 0.9]';
s = [2.7, 1.48, 0.819, 0.458]';
x0 = [1, 1, 1, 2]';
rj = @(j, x) sum(arrayfun(@(k) x(k)*exp(-x(n+k)*t(j)), 1:n)) - s(j);
r = @(x) arrayfun(@(j) rj(j,x),1:2*n)';
Jj1 = @(j, x) arrayfun(@(i) exp(-x(n+i)*t(j)),1:n);
Jj2 = @(j, x) arrayfun(@(i) -t(j)*x(i-n)*exp(-x(i)*t(j)),n+1:2*n);
Jj = @(j, x) [Jj1(j, x)'; Jj2(j, x)'];
J = @(x) cell2mat(arrayfun(@(j) Jj(j, x), 1:2*n, 'UniformOutput', false))';

TOL = 1e-10;
MaxIter = 1000;

x = gaussNewton(x0, J, r, TOL, MaxIter);
k = length(x);
sol = x(:,end);

f = @(t) sum(arrayfun(@(k) sol(k)*exp(-sol(n+k)*t), 1:n));
plot(linspace(0,1,100), arrayfun(f, linspace(0, 1, 100)), 'k-', t, s, 'ro');
