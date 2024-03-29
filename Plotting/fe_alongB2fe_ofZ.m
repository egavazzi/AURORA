function fe_alt = fe_alongB2fe_ofZ(h_atm,feZE,iZ)
%  fe_alongB2fe_ofZ - extract energy-pitch-angle phase-space-density 
%  plot_fezv_2Dvtheta extracts the pitch-angle energy electron
%  phase-space-density  at a selected altitude. 
%
% Calling:
%  output = fe_alongB2fe_ofZ(h_atm,feZTE,BeamW,iZ,movieOut)
% Input:
%  h_atm      - altitudes (km), double array [n_z x 1]
%  fe_ztE     - electron number-flux, double array [n_z*n_beams,n_t,n_E]
%  iZ         - index in altitude, scalar integer, iZ should be in
%               [1,...,n_z]

%   Copyright � 20190506 Bj�rn Gustavsson, <bjorn.gustavsson@uit.no>
%   This is free software, licensed under GNU GPL version 2 or later




nZ = numel(h_atm);
%% Extracting electron flux at 4 altitudes


fe_alt = feZE(iZ:nZ:end,:);
