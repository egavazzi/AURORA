function Ie = Ie_flickeringGauss_flux(t,E,I0,E0,dE,Ie_LET,E_LET,w,z0,Emax,Iflicker,Ef,dEf,mu,mu0)
% IE_FLICKERINGGAUSS_FLUX - Gaussian with LET harmonically modulated
% electron spectra. IE_FLICKERINGGAUSS_FLUX gives electron spectra
% with Gaussian number-flux with peak flux at E0 and a width in
% energy of dE with a peak electron number flux of I0. The function
% handles an energy dependent low-energy tail, if given as input
% parameter. A Gaussian flux with peak flux Iflicker at energy Ef
% with width dEf is modulated at w (angular frequency) with a
% sin(w*t/2)^2 wave at all energies at a distance of z0 above,
% allowing for time-dispersion due to collisionless transport. This
% spectra is suitable to model flickering electron precipitation at
% high enegries. This to allow for comparisons with low-energy
% FAB-flickering with harmonic and on-off modulated flux and their
% corresponding ionospheric responses. 
% 
% Calling:
%   Ie = Ie_flickeringGauss_flux(t,E,I0,E0,dE,Ie_LET,E_LET,w,z0,Emax,Iflick,Ef,dEf,mu,mu0)
% Input:
%  t      - time (s), double array [1 x nT]
%  E      - energy (eV), scalar double
%  I0     - peak electron number flux of Gaussian (#e/m^2/s) in an
%           energy bin
%  E0     - centre-energy of Gaussian (eV), double scalar
%  dE     - energy-width of Gaussian (eV),
%  Ie_LET - differential electron flux (#e^-/m^2/s/eV) of
%           low-energy tail, double array [1 x nE]
%  E_LET  - energy (eV) grid of LET, double array [1 x nE]
%  w      - angular frequency of modulation (radians/s), double
%           scalar
%  z0     - distance (m) above "top-of-ionosphere", double scalar
%  Emax   -
%  Iflick - peak flickering electron number flux of Gaussian
%           (#e/m^2/s) in an energy bin
%  Ef     - centre-energy of flickering Gaussian (eV), double scalar
%  dEf    - energy-width of flickering Gaussian (eV),
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

% Handy ruler for nargin-counting
% IE_FLICKERINGGAUSS_FLUX -           1,2, 3, 4, 5,     6,    7,8, 9,  10,      11,12, 13,14, 15


if nargin < 11 || isempty(Iflicker)
  Iflicker = 0;
end
if nargin < 12 || isempty(Ef)
  Ef = E0;
end
if nargin < 13 || isempty(dEf)
  dEf = dE;
end
if nargin < 14 || isempty(mu)
  mu = 1;
end
if nargin < 15 || isempty(mu0)
  mu0 = 1;
end

t_shift0 = z0/(mu0*v_of_E(Emax));
t_shift  = z0./(mu*v_of_E(E));
t_shifted = t-(t_shift-t_shift0);

Ie = I0*exp(-(E-E0).^2/dE^2) + ...
     Iflicker*exp(-(E-Ef).^2/dEf^2)*sin(w/2*t_shifted).^2.*(t_shifted>0);
[~,iE] = min(abs(E-E_LET));

Ie = max(Ie, Ie_LET(iE));
