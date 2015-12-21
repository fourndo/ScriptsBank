function derivTest(f, df, n)
%DERIVTEST Performs derivative test on simulator operators f and df.
%   f should be a function handle of argument in one variable x
%   of size n x 1 (row vector length n). It should perform the
%   forward simulation (or generically speaking the matrix operation on x)
%   and output a size(x) vector d (data) (length n row vector).
%   df should be the derivative operator of f, and should therefore be size
%   n x n. ie take a vector of length n and output the derivative length n.
%   
%   The Derivative test is a crude test of the accuracy of matrix operators
%   and boundary conditions. It relies on the Taylor Expansion eq 1:
%   1) f(x+a) = f(x) + af'(x) + O(a^2)
%   2) f(x+a) - f(x) = O(h)
%   3) f(x+a) - f(x) - af'(x) = O(a^2)
%
%   As h is decreased by orders of magnitude equation 2) should drop of as
%   O(h) and equation 3) should drop off as O(h^2).

% Random model
x = randn(n, 1);
% Random perturbation f(x+a) where a = h*v
v = randn(n, 1);

fx = f(x);
dfdx = df(x);

N = 10;
% decrease h by an order of magnitude each iteration
% Take norm of left hand side of eqn 2) and 3), take slope, should act
% as right hand side of 2) and 3).

for jj = 1 : 10
    h = 10^(-jj);

    err2 = norm(f(x + h*v) - fx); %#ok<*AGROW>
    err3 = norm(f(x + h*v) - fx - h * dfdx * v);

    fprintf('%3.2e    %3.2e    %3.2e\n',h , err2, err3)
    
end
end
