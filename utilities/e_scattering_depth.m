function tau = e_scattering_depth(varargin)
% E_SCATTERING_DEPTH - electron scattering depth
%   
% Calling:
%  tau = e_scattering_depth(z,nA,xsA,nB,xsB,...,nX,xsX)
% Input:
%  nA  - number-density of species A, double array [nA x 1]
%  xsA - collision cross-section array, double array [nLevelsA x nE],
%        all collision cross-sections for electrons with
%        atoms/moelcules of species A, ellastic, inellastic,
%        ionization etc.
%  nB  - number-density of next species, double array [nA x 1]
%  xsB - collision-cross-sections for electrons with
%        atoms/molecules of species B, double array [nLevelsB x nE]
%  :
%  :     and so on...
% Output:
%  l_meanfree - mean-free-path, double array [nA x nE]
% 
% Units should match to give sensible dimension of mean-free-path,
% i.e. m^-3 and m^2 for number-density and cross-section giving
% l_meanfree the unit of m.
% 
% The function will calculate the mean-free-path at all levels of
% the neutral density, typically altitude-profiles at all energies
% giving a 2-D array of mean-free-paths. 

% Copyright © 20190502 B. Gustavsson, <bjorn.gustavsson@uit.no>
%   This is free software, GPL version 3 or later applies

dSC = 0;

for i1 = 2:2:nargin
  
  n_i = varargin{i1};   % Neutral density
  xs_i = varargin{i1+1};% Array with colission cross sections
  
  dSC = dSC + n_i*xs_i(1,:);
  % add the inelastic colissions to the total
  for j2 = 2:size(xs_i,1)-1,
    
    dSC = dSC + n_i*xs_i(j2,:);
    
  end
  
end

z = varargin{1};
tau = -flipud(cumtrapz(flipud(z),flipud(dSC)));
