function Iei = Ie_tmuE_par(t,E,mu,Ie,ti,Ei,mui,Eref,dt_of_mu,Mu_at_source,dt00,dmu,dE,method)
% Ie_tmuE_par - energy-pitchangle-time-dispersion electron-fluxes reinterpolation
% Ie_tmuE_par reinterpolates electron-fluxes taking into account
% energy and pitch-angle time-of-fligth dispersion from one
% altitude to an altitude-of-interest. The time-shift for electrons
% with different pitch-angle-cosines at the target-altitude have to
% be given, as will the pitch-angle-mapping between the source and
% target altitudes.
%
% Calling:
%  Iei = Ie_tmuE_par(t,E,mu,Ie,ti,Ei,mui,Eref,dt_of_mu,Mu_at_source,dt00,dmu,dE,method)
% Input:
%  t      - time (s) of electron spectra at source altitude, double
%           array [1 x nt]
%  E      - energy (eV) of electron spectra at source altitude,
%           double array [nE x 1]
%  mu     - pitch-angle-cosine at source altitude, double array 
%           [1 x n_mu], typically cosd(0:89)
%  Ie     - Electron flux array [e^-/eV/m^2/s], double array [nE x nt]
%  ti     - time-array (s) to evaluate electron-fluxes at, double
%           array [1 x n_ti], typically at finer resolution than t
%  Ei     - energy-array (eV) to evaluate electron-fluxes at, double
%           array [n_Ei x 1], typically at finer resolution than E
%  mui    - Pitch-angle-cosine at the target altitude to calculate
%           the electron-fluxes for 
%  Eref   - Reference-energy (eV) that the time-shifts in dt_of_mu
%           are calculated for.
%  dt_of_mu - Time-of-flights (s) between the source-altitude and the
%             target altitude for electrons with energy Eref for
%             pitch-angle-cosines mu at the target-altitude, double
%             array [1 x n_mu]
%  Mu_at_source - pitch-angle-cosine-array with pitch-angle-cosines
%                 at the source altitude that end up with
%                 pitch-angle-cosines mu at the target altitude,
%                 double array [1 x n_mu]
%  dt00   - base-time-shift (s), double scalar, used to shift the
%           shortest time-of-flight to zero. Double scalar, set to
%           zero to disable
%  dmu    - pitch-angle-cosine-width of electron-fluxes at
%           source-altitude 
%  dE     - energy-element width (eV)
%  method - interpolation-method, string [{'linear'} | 'nearest']
% Output:
%  Iei - electron-fluxes (e^-/m^2/s/dE) at target altitude, double
%        array [n_Ei x n_ti]
%  
% 
% Example:
%  load('Ie.mat','Ec','tC','Ie_Bpar')
%  Eref = max(Ec);
%  dE = 10;
%  [B,lB,lon,lat,alt,x_B,y_B,z_B] = make_BofL([],[],[],h_atm(end),z_0);
%  for i1 = 1:90,
%    [t_delay,mu_source] = dt_of_z_mirror(v_of_E(Eref),lB(1),lB(end),cos((i1-1)*pi/180),B,lB);
%    dt0(i1) = t_delay;                                                                     
%    MU_at_source(i1) = mu_source;                                                          
%  end
%  dt00 = min(dt0);
%  for imui = numel(mui):-1:1,
%    for iE = numel(E):-1:1,    
%     Ie_test(iE,imui,:) = Ie_tmuE_par(tC,Ec,cos((0:89)*pi/180),...
%                                      Ie_Bpar,...
%                                      ti,E(iE),mui(imui),...
%                                      Eref,dt0,MU_at_source,...
%                                      dt00,cos(9*pi/180),10,'linear'); 
%    end
%  end
%  subplot(4,1,1)                                                              
%  pcolor(ti,E,log10(squeeze(Ie_test(:,1,:)))),
%  shading flat,caxis([-6 0]+max(caxis))
%  subplot(4,1,2)
%  pcolor(ti,acos(mui)*180/pi,log10(squeeze(Ie_test(80,:,:)))),shading flat
%  subplot(4,1,3)                                                              
%  pcolor(ti,acos(mui)*180/pi,log10(squeeze(Ie_test(40,:,:)))),shading flat
%  subplot(4,1,4)                                                          
%  pcolor(ti,acos(mui)*180/pi,log10(squeeze(Ie_test(20,:,:)))),shading flat

%   Copyright © 2020 Bjorn Gustavsson, <bjorn.gustavsson@uit.no>
%   This is free software, licensed under GNU GPL version 2 or later

% input-argument-ruler    1 2  3  4  5  6   7    8        9           10   11   12 13     14     15

% Set default fallbacks
if nargin < 13 || isempty(dE)
  dE = 1;
end
if nargin < 14 || isempty(method)
  method = 'linear';
end

% calculate the pitch-angle-cosine at the source-altitude for
% pitch-angle-cosine mui at the target altitude.
mu_at_source = interp1(mu,Mu_at_source,mui);

for iE = 1:numel(Ei)
  % Energy-scaled time-of-flight between source and target
  % altitudes for electrons with target pitch-angle-cosine mui,
  % shifted with dt00.
  t_thetashift = interp1(mu,dt_of_mu,mui)*v_of_E(Eref)/v_of_E(Ei(iE)) - dt00;
  
  % time-shifted electron-flux re-interpolation with pitch-angle-width
  Iei(iE,:) = interp2(t,...
                      E,...
                      Ie,...
                      ti-t_thetashift,...
                      Ei(iE),...
                      method,0)*exp(-abs(acos(mu_at_source)).^2/acos(dmu)^2).*dE;
end
Iei(~isfinite(Iei(:))) = 0;
