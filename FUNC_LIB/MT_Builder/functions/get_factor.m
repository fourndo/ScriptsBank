function [a]=get_factor(n,dx,L)
% Function computes the coefficiant for the expanding vector. 
% [a]=get_factor(n,dx,L);
% The length (L) of the vector can be written as:
% L=dx + dx*a + dx*a^2 + ... + dx*a^n, 
% where dx is the initial cell size
%       a is the expansion factor
%       n is the number of cells to reach L
%
% The function solves for the expansion coefficient by computing the roots
% of the polynomial. The only real and positive root is returned
% Written by: D. Fournier
% Last Update: 2013-02-01

polynom=ones(1,n)*dx;
polynom(end)=polynom(end)-L;
getroots=roots(polynom);
a=(real(getroots(end)));