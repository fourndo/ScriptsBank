function Derivative_test(C_p,dCdm_p,L,nface,mcell)
% Derivative test for primary field operator using random vectors
% From the derivative of the forward operation, the fuunction computes:
% f(x+h) = f(x) + f'(x) * h

a = randn(nface, 1);
phi = randn(mcell, 1);
u = [a;phi];
m = randn(mcell, 1);
q = randn(mcell+nface, 1);
v = randn(mcell, 1);

fprintf('**Derivative Test**\n')
fprintf('\t1/h\t\t|f(x)-f(x+h)|\t|f(x+h) - f(x)- dfdx*h|\n')

for ii=1:10
h = 10^(-ii);

m2=( m + v*h );


Oh = norm(C_p(m2,u,q,L) - C_p(m,u,q,L));
Oh2 = norm(C_p(m2,u,q,L) - C_p(m,u,q,L) - dCdm_p(a,phi,m)*h*v);

fprintf('%3.2e    %3.2e \t\t\t %3.2e\n',h,Oh,Oh2)
end
fprintf('**end of test**\n\n')
