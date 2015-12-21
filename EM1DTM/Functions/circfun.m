function [x, y] = circfun(xc, yc, r, n)
theta = linspace(-pi, pi, n+1);
x = r*cos(theta(:))+xc;
y = r*sin(theta(:))+yc;

x = x(1:end-1);
y = y(1:end-1);
end