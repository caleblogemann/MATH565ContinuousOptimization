% Discretization of the 5-point stencil for the lexicographically ordered
% Laplace's equation on the 2-D unit square.
%
% Discretization of
%    u_xx + u_yy = 0  on  alpha<=x<=beta and gamma<=y<=delta.
%    u(alpha,y) = f1
%    u(beta,y)  = f2
%    u(x,gamma) = f3
%    u(x,delta) = f4
% Using n number of points with equal step size h1 in both x and y.

function A = CreateA(n1)

    % Parameters which can be varied
    alpha=0; beta=1; gamma=0; delta=1;
    f1=0; f2=0; f3=0; f4=0;

    % 'A' Matrix Composition
    h1=(beta-alpha)/(n1+1);

    T = sparse(n1,n1);
    I = sparse(n1,n1);
    e1=ones(n1,1);
    e2=zeros(n1,1);
    T=spdiags([e1 -4*e1 e1], [-1 0 1], n1, n1);
    I=spdiags([e1],[0],n1,n1);

    A = sparse(n1*n1,n1*n1);
    for i=1:n1
        
        ma=(i-1)*n1+1;
        mb=i*n1;
        
        A(ma:mb,ma:mb)=T;
        
    end

    for i=1:(n1-1)
        
        ma = (i-1)*n1+1;
        mb = i*n1;
        
        mc = i*n1+1;
        md = (i+1)*n1;
        
        A(ma:mb,mc:md)=I;
        A(mc:md,ma:mb)=I;
        
    end

    A = -A;
end