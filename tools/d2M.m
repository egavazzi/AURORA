function D2M = d2M(z)
% D2M - second-order differentiation matrix for non-uniform grid
%   D2M calculates the second-order differentiation-matrix for data
%   on a non-uniform one-dimensional grid. This is the properly
%   weighted differential operator centred at the grid-points. The
%   first and last rows of the differentiation matrix is set to
%   zero, so the first and last element of the second derivatives
%   are returned as zeroes.
% 
% Calling:
%   D2M = d2M(z)
% Input:
%   z   - grid coordinate vector, double array [n_z x 1] or [1 x n_z]
% Output:
%   D2M - differentiation matrix, double array [n_z x n_z]
% Example:
%  x = [1.9;2;2.2;2.5]; % non-uniform grid, advanced stuff!
%  D2M = d2M(x)
%  y = x.^2;
%  d2ydx2 = D2M*y
%  % analytical solution:
%  % dydx = 2*x
%  % d2ydx2 = 2
%  y = x.^3;
%  d2ydx2 = D2M*y
%  % Analytical solution
%  % dydx = 3*x^2
%  % d2ydx2 = 2*3*x
%  % Comparison:
%  [d2ydx2(2:end-1)./(2*3*x(2:end-1))]
%  % First order accurate in gradient-of-grid-spacing, see below
%
% No argument checks or error-controls, if you want that you have
% to pay me good money.

% Copyright © 20190506 B. Gustavsson, <bjorn.gustavsson@uit.no>
%   This is free software, GPL version 3 or later applies

dzd = z(2:end-1) - z(1:end-2);
dzu = z(3:end)   - z(2:end-1);


dsup  = [2./(dzd(:).*(dzu(:)+dzd(:)));0];
dMain = [0;-2./(dzd(:).*dzu(:));0];
dsub  = [0;2./((dzd(:)+dzu(:)).*dzu(:))];

D2M = diag(dMain) + diag(dsub,1) + diag(dsup,-1);

%%  Derivation using the matlab symbolic toolbox
% syms dxu dxd    
% syms f fu fd    
% syms dfdx d2fdx2
% qwe = solve(fu == f + dxu*dfdx + (dxu)^2/2*d2fdx2,fd == f - dxd*dfdx + (dxd)^2/2*d2fdx2,dfdx,d2fdx2);
% qwe.d2fdx2 = -(2*(dxd*f + dxu*f - dxu*fd - dxd*fu))/(dxd*dxu*(dxd + dxu))
% pretty(simplify(expand(qwe.d2fdx2)))
%
%    2 (dxd f + dxu f - dxu fd - dxd fu)
%  - -----------------------------------
%            dxd dxu (dxd + dxu)
%Manual reordering
%
%    2 dxd f + 2 dxu f - 2 dxu fd - 2 dxd fu
%  - ---------------------------------------
%              dxd dxu (dxd + dxu)
%
%     - 2 dxu fd + 2 dxd f + 2 dxu f - 2 dxd fu
%  - ---------------------------------------
%              dxd dxu (dxd + dxu)
%
%      2 dxu fd - 2 (dxd + dxu) f + 2 dxd fu
%    ---------------------------------------
%              dxd dxu (dxd + dxu)
%
%      2 dxu fd             2 (dxd + dxu) f       2 dxd fu
%    ------------------- - ------------------ + ---------------
%    dxd dxu (dxd + dxu)   dxd dxu (dxd + dxu)  dxd dxu (dxd + dxu)
%
%      2 fd             2  f       2  fu
%    --------------- - -------- + ---------------
%    dxd (dxd + dxu)   dxd dxu    dxu (dxd + dxu)
% 
% % Error-order-estimate:
% syms d3fdx3
% asd = solve(fu == f + dxu*dfdx + (dxu)^2/2*d2fdx2 + (dxu)^3/6*d3fdx3,...
%             fd == f - dxd*dfdx + (dxd)^2/2*d2fdx2 - (dxd)^3/6*d3fdx3,dfdx,d2fdx2);
%
% Leading error term:
%
% pretty(simplify(expand(asd.d2fdx2 - qwe.d2fdx2)))
%
%  d3fdx3 (dxd - dxu)
%  ------------------
%          3
