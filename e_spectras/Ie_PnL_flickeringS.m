function Ie = Ie_PnL_flickeringS(t,E,IeEtot,Ein,dEin,w,z0,mu,mu0)
% IE_PNL_FLICKERINGS - On-off modulated P-n-L flickering spectra
% for electrons. Ie_PnL_flickeringH produces electron spectra with
% constant number-flux between the lowest and highest energy of Ein
% with an integrated energy flux of IeEtot. The flux is modulated
% at w (angular frequency) with a square(w*t-pi/2)^2) wave at all
% energies at a distance of z0 above, allowing for time-dispersion
% due to collisionless transport. This Peticolas and Lummerzheim
% [2000] spectra is modified to have its modulation in phase with
% Ie_PnL_flickeringH making it suitable to model low-energy-tail  
% field-aligned flickering precipitation with a harmonically
% modulated intensity. This to allow for comparisons of harmonic
% and on-off modulated flux and their corresponding ionospheric
% responses. 
% 
% Calling:
%  Ie = Ie_PnL_flickeringS(t,E,IeEtot,Ein,dEin,w,z0,mu,mu0)
% Input:
%  t      - time (s), double array [1 x nT]
%  E      - energy (eV), scalar double
%  IeEtot - total energy flux (W/m^2) of the FAB
%  Ein    - energy grid of the FAB (eV), double array [1 x nE]
%  dEin   - size of energy bins (eV), same size as Ein
%  w      - angular frequency of modulation (radians/s), double
%           scalar
%  z0     - distance (m) above "top-of-ionosphere", double scalar
%  mu     - pitch-angle average of precipitation (radians), double
%           scalar used to calculate time-lag of precipitation
%           relative to electrons at energy max(Ein) with
%           pitch-angle mu0
%  mu0    - reference pitch-angle (radians), scalar double
% Output:
%  Ie - differential electron flux (#/m^2/s), double array [1 x nT]
% 
% NOTE flux is per energy bin!

% Copyright © B Gustavsson 20180512, bjorn.gustavsson@uit.no
% This is free software, licensed under GNU GPL version 2 or later

 



q_e     = 1.602176620898e-19;

t_shift0 = z0/(mu0*v_of_E(Ein(end)));
t_shift  = z0./(mu*v_of_E(E));

t_shifted = t-(t_shift-t_shift0);

IePnL = IeEtot/q_e/sum((Ein+dEin/2).*dEin).*dEin;
[~,iE] = min(abs(E-Ein));

Ie = IePnL(iE)*(1+square(w*t_shifted-pi/2))/2.*(t_shifted>0)*double(Ein(1)<=E)*double(E<=Ein(end));