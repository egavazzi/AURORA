function cblh = plot_fezv_2DvCthetaE(h_atm,E,dE,feZE,BeamW,mu_lims,iZ,cxin)
%  plot_fezv_2DvCthetaE - energy-pitch-angle phase-space-density plot
%  plot_fezv_2Dvtheta plots the pitch-angle energy electron
%  phase-space-density plots at a selected altitude. The
%  phase-space-density is plotted in log-scale in energy with pcolor.
%
% Calling:
%  output = plot_fezv_2Dvtheta(h_atm,E,dE,feZTE,BeamW,mu_lims,iZ,movieOut)
% Input:
%  h_atm      - altitudes (km), double array [n_z x 1]
%  E          - energies (eV), double array [1 x n_E]
%  fe_ztE     - electron number-flux, double array [n_z*n_beams,n_t,n_E]
%  dE         - energy bin sizes (eV), double array [1 x n_E]
%  BeamSA     - solid angle sizes of the beams, [n_beams x 1]
%  mu_lims    - beam-boundaries, (cos(theta_b)), [n_beams x 1] or
%               [1 x n_beams]
%  iZ         - index in altitude, scalar integer, iZ should be in
%               [1,...,n_z]

%   Copyright © 20190506 Björn Gustavsson, <bjorn.gustavsson@uit.no>
%   This is free software, licensed under GNU GPL version 2 or later


nZ = numel(h_atm);
%% Extracting electron flux at 4 altitudes


fe_alt = feZE(iZ:nZ:end,:);

fe_max = log10(max(fe_alt(:)));
%% Graphics settings:
ftsz = 14;
xtxtpos = [0,1.3]; % Altitude-label above change 1.3 to -1.3 for below
%suptpos = [1.1,1.4];
cblpos = [-0.7,-0.3];

%% Animationing

pcolor((E(1:size(feZE,2))),...
       acos(mu_lims)*180/pi,...
       max(fe_max-10,log10(fe_alt([1:end,end],:))));
shading flat
set(gca,'xscale','log')
ylabel('Pitch-angle (deg)','fontsize',ftsz);
xlabel('Energy (eV)','fontsize',ftsz);

alt_str = sprintf('z: %2.1f (km)',h_atm(iZ));
title(alt_str,'fontsize',ftsz);
clbh = colorbar_labeled('(e^{-1}s^3/m^6)','log','fontsize',ftsz-2);

