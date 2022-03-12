function Ie = Ie_E_with_LET(E0,Q,E,LEToff)
% IE_E_WITH_LET - electron spectra with low energy tail (energetic e)
%   
% IE_E_WITH_LET gives the flux from a Maxwellian spectra with a low
% energy tail - implemetation of Meier/Strickland/Hecht/Christensen
% JGR 1989 (pages 13541-13552)
% 
% Calling:
%  Ie = Ie_E_with_LET(E0,Q,E)
% Input:
%  E0 - characteristic energy (eV)
%  Q - energy flux (eV/cm^2/s)
%  E - energy grid
% Output:
%  Ie - differential electron flux (#/eV/m^2/s/ster)
% 
% NOTE flux is per eV!!!

% Copyright © B Gustavsson 20070502, bjorn.gustavsson@uit.no
% This is free software, licensed under GNU GPL version 2 or later


% Parameter for LET amplitude
b = (0.8*E0/1e3).*(E0<500) + (0.1*E0/1e3+.35).*(E0>=500);

% Maxwellian spectra
Ie = Q/(2*pi*E0^3).*E.*exp(-E/E0);


if nargin > 3 && LEToff
else
  % Max of Maxwellian - to scale LET amplitude
  Ie_max = max(Ie);
  
  Ie = Ie + 0.4*Ie_max*(E0./E).*exp(-E/1e3/b);
  
end
