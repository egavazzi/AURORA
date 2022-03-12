function [f_e,v_E] = Ie_pad2fe(Ie_E,E,BeamSA)
% IE_PAD2FE - Electron flux to phase-space density conversion
% Ie_pad2fe converst multi-stream electron fluxes to phase-space
% densities.
%   
% Calling:
%  [f_e,v_E] = Ie_pad2fe(Ie_E,E,BeamSA)
% Input:
%  Ie_E   - Multi-stream electron fluxes, (e^-1/m^2/s/DE), double
%           array, either [(n_z*n_mu) x n_E] or [(n_z*n_mu) x n_t x n_E] 
%           where n_z is the number of altitude levels, n__mu is
%           the number of pitch-angle-streams, and n_E is the
%           number of energy-bins. Ie_E can be a time-resolved
%           electron flux modeling the time-variatino of
%           electron-fluxes, or a steady-state electron-flux (or a
%           single-time-step snap-shot). 
%  E      - Energy (eV), double array [1 x n_E]
%  BeamSA - solid-angles of pitch-angle-streams (ster), double
%           array [n_mu x 1] or [1 x n_mu]
% Output:
%  f_e - phase-space-density (e^-1 s^3/m^6), double array with the
%        same dimensions as Ie_E. f_e are discretised in the same
%        way as Ie_E, that is in n_mu pitch-angle-streams, with nE
%        steps in energy
%  v_E - velocity (m/s) at centre of energy bin, double array [1 x n_E]

%   Copyright © 2020200115 Björn Gustavsson, <bjorn.gustavsson@uit.no>
%   This is free software, licensed under GNU GPL version 2 or later

n_mu = numel(BeamSA);
sz_Ie = size(Ie_E);
nE = sz_Ie(end);
nZ = sz_Ie(1)/n_mu;

dE = diff(E);
dE(end+1) = dE(end);

if numel(sz_Ie) == 3 % the Ie_ztE case
  v_E(1,1,:) = v_of_E(E+dE/2);
  dv_E(1,1,:) = v_of_E(E+dE) -v_of_E(E);
else % the Ie_zE case
  v_E = v_of_E(E+dE/2);
  dv_E = v_of_E(E+dE) -v_of_E(E);
end

%% Conversion of flux to phase-space density
% This takes us to phase-space density in units of /m^3/d^3v
% that is per discrete velocity cell (since the input electron flux
% it supposed to be in [e^-/m^2/s/dE] (electrons per square metre
% per second in each energy bin))
f_e = Ie_E./repmat(v_E,[sz_Ie(1:end-1),1]);

% Continuing converting to proper phase-space-density units
% (/m^3/(m/s)^3)
f_e = f_e./repmat(v_E.^2.*dv_E,[sz_Ie(1:end-1),1]);

for i_mu = 1:n_mu,
  f_e(nZ*(i_mu-1) + (1:nZ),:) = f_e(nZ*(i_mu-1) + (1:nZ),:)/BeamSA(i_mu);
end
