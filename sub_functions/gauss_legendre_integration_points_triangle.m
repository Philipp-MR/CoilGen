function [u,v,ck] = gauss_legendre_integration_points_triangle(n)
% This file aim at calculating the coordinates needed for a Gauss Legendre
% integration method on a triangle and the associated weighting factors
% n : order of the foreseen integration
% u,v : coordinate of the point
% ck : weighting coefficients

[eta,w] = calc_weights_gauss(n);

k=1;
for i=1:size(eta,1)
    for j=1:size(eta,1)
        u(k,1) = (1+eta(i))/2;
        v(k,1) = (1-eta(i))*(1+eta(j))/4;
        ck(k,1) = ((1-eta(i))/8)*w(i)*w(j);
        k=k+1;
    end
end


function [g_abscissa, g_weights] = calc_weights_gauss(n)
% Generates the abscissa and weights for a Gauss-Legendre quadrature.
% Reference:  Numerical Recipes in Fortran 77, Cornell press.
g_abscissa = zeros(n,1);                                           % Preallocations.
g_weights = g_abscissa;
m = (n+1)/2;
for ii=1:m
    z = cos(pi*(ii-.25)/(n+.5));                        % Initial estimate.
    z1 = z+1;
while abs(z-z1)>eps
    p1 = 1;
    p2 = 0;
    for jj = 1:n
        p3 = p2;
        p2 = p1;
        p1 = ((2*jj-1)*z*p2-(jj-1)*p3)/jj;       % The Legendre polynomial.
    end
    pp = n*(z*p1-p2)/(z^2-1);                        % The L.P. derivative.
    z1 = z;
    z = z1-p1/pp;
end
    g_abscissa(ii) = -z;                                   % Build up the abscissas.
    g_abscissa(n+1-ii) = z;
    g_weights(ii) = 2/((1-z^2)*(pp^2));                     % Build up the weights.
    g_weights(n+1-ii) = g_weights(ii);
end
end

end