function [x,w]=lgwt(N,a,b)
% [x0,w] = LGWT(N,a,b) - Computes the Legendre-Gauss nodes and weights on an interval [a,b]
%   with truncation order N.
%
%  The definite integral in interval [a,b] for a continuous function f(x) approximated by a
%  polynomial of order N is given by sum(f(x0).*w);
%
% MIT License, Copyright (c) 2016 Pazus <https://github.com/Pazus/Legendre-Gauss-Quadrature>
% Based on Greg von Winckel (2004) <https://www.mathworks.com/matlabcentral/fileexchange/...
%   4540-legendre-gauss-quadrature-weights-and-nodes>

if nargin == 1
    a = -1;
    b =  1;
end
x = zeros(1, N);
w = zeros(1, N);
m = (N+1)/2;
h = b-a;

for ii=1:m
    z = cos(pi*(ii-.25)/(N+.5)); % Initial estimate.
    z1 = z+1;
    while abs(z-z1)>eps
        p1 = 1;
        p2 = 0;
        for jj = 1:N
            p3 = p2;
            p2 = p1;
            p1 = ((2*jj-1)*z*p2-(jj-1)*p3)/jj; % The Legendre polynomial.
        end
        pp = N*(z*p1-p2)/(z^2-1); % The L.P. derivative.
        z1 = z;
        z = z1-p1/pp;
    end
    x(ii) = z; % Build up the abscissas.
    x(N+1-ii) = -z;
    w(ii) = h/((1-z^2)*(pp^2)); % Build up the weights.
    w(N+1-ii) = w(ii);
end

if a ~= -1 || b ~= 1
    x = (x+1)*(h/2) + a;
end
