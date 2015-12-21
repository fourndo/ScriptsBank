function F = forwardOp(Obs1, Obs2, Obs3, xc1, xc2, xc3, h)
% FORWARDOP calculate the forward operator. Note depth dimension is X1

% Get number of cells in each direction
n1 = numel(xc1);
n2 = numel(xc2);
n3 = numel(xc3);

% Get numbor of total cells and observations
nnn = n1 * n2 * n3;
nobs = numel(Obs1);

% Helper funcs
e = @ (n) ones(1, n);
% this krons from 1D -> 3D
kron3 = @(a, b, c)  kron(a, kron(b, c));

% Kron observations to size of Forward Operator
kronobs = @(obs) kron( e(nnn), obs(:) );

% Kron dimensions to size of forward OP and subtract obs distance
H = kron( e(nobs)', h');
X1 = kron( e(nobs)', kron3( e(n3), e(n2), xc1')) - kronobs(Obs1);
X2 = kron( e(nobs)', kron3( e(n3), xc2', e(n1))) - kronobs(Obs2);
X3 = kron( e(nobs)', kron3( xc3', e(n2), e(n1))) - kronobs(Obs3);
R = (X1.^2 + X2.^2 + X3.^2).^(3/2);

F = X1 .* H .* 1./R;

end