function Le = dEdt_ee(E,ne,Te)
% dEdt_ee - suprathermal electron energy loss rate in e-e-collisions 
%  This function is the electron energy loss rate due to
%  electron-electron interaction. This is the analytic form given
%  for the energy-transfer rate from photoelectrons to thermal
%  electrons, given by Schwartz and Nisbeth (1971). The expression
%  fits the classical formulation of Itakawa and Aono (1966) at low
%  energies and gives a  smooth transition to fit the quantum
%  mechanical equation of Schunk and Hays (1971) at higher
%  energies. 
% 
% Calling:
%  Le = dEdt_ee(E,ne,Te)
% Input:
%  E  - energy levels [1 x nE] (eV), 
%  ne - ambient electron concentration  [nZ x 1] (/m^3),
%  Te - electron temperature [nZ x 1] (K).
% Output:
%  Le - electron energy loss-rate [eV/s], [nZ x nE]
% 
% The paper by Schwartz and Nesbith use electron density in cm^-3,
% here the constant is rescaled to use m^-3 instead.
% 
% References:
%  * Swartz, W. E., J. S. Nisbet, and A. E. S. Green (1971), Analytic
%  expression for the energy‐transfer rate from photoelectrons to
%  thermal‐electrons, J. Geophys. Res., 76(34), 8425-8426,
%  doi: 10.1029/JA076i034p08425.
%  * Itikawa, Y., and O. Aono (1966), Energy change of a charged particle
%  moving in a plasma, Phys. Fluids, 9, 1259-1261.
%  * Schunk, R. W., and P. B. Hays (1971), Photoelectron-energy
%  losses to thermal electrons, Planet. Space Sci., 19, 113-117.

% Copyright ,A) B. gustavsson 20180527
%  This is free software, licensed under GNU GPL version 2 or later

kB	= 1.380662e-23;		    % Boltzmann constant [J/K]
q_e     = 1.6021773e-19;            % elementary charge [C]

if ( length(E) > 1 )
  
  sz_E = size(E);
  sz_ne = size(ne);

  ne = repmat(ne,sz_E);% *e_dim;
  Te = repmat(Te,sz_E);% *e_dim;
  E = repmat(E,sz_ne); % alt_dim*E;
  Ee = kB/q_e*Te;
  
  Le = real(3.0271e-10*ne.^.97.*((E-Ee)./(E-.53*Ee)).^2.36./E.^.44);
  
  Le = Le.*double(E>=Ee);
else
  
  Ee = 8.618e-5*Te;
  
  Le = real(3.0271e-10*ne.^.97.*((E-Ee)./(E-.53*Ee)).^2.36./E.^.44);
  Le(E<Ee) = 0;
  
end
