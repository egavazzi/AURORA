function output = extract_IeztE_pitchangleflux(t,h_atm,dE,IeZTE,BeamW,iZ)
% animate_IeztE_pitchangleflux - animation of time-varying electron energy flux
%  animate_IeztE_pitchangleflux animates the time-varying electron
%  fluxes varying with energy - pitch-angle at four altitudes
%  stepping in time. The electron-fluxes are plotted in log-scale
%  with polarPcolor in four subplots.
%
% Calling:
%  output = animate_IeztE_pitchangleflux(t,h_atm,E,dE,IeZTE,BeamW,mu_lims,iZ,movieOut)
% Input:
%  t          - time-scale (s), double array [1 x n_t]
%  h_atm      - altitudes (km), double array [n_z x 1]
%  E          - energies (eV), double array [1 x n_E]
%  Ie_ztE     - electron number-flux, double array [n_z*n_beams,n_t,n_E]
%  dE         - energy bin sizes (eV), double array [1 x n_E]
%  BeamSA     - solid angle sizes of the beams, [n_beams x 1]
%  mu_lims    - beam-boundaries, (cos(theta_b)), [n_beams x 1] or
%               [1 x n_beams]
%  movieOut   - Optional argument, if MovieName in matlab-supported
%               video-format then animate_IeEztE_3DtEofz will attempt
%               to write the displayed sequence to a video-file
%               with that filename using the VideoWriter
%               functionality. See VideoWriter for details. If
%               numerical value that is converted to TRUE, a matlab
%               movie will be returned.
% Output:
%  M     - Matlab-movie with the displayed sequence.


%   Copyright © 20190506 Björn Gustavsson, <bjorn.gustavsson@uit.no>
%   This is free software, licensed under GNU GPL version 2 or later



nZ = numel(h_atm);
%% Extracting electron flux at 4 altitudes

for iz = numel(iZ):-1:1
  Ie_alts{iz} = IeZTE(iZ(iz):nZ:end,:,:);
end

%% Graphics settings:

%% Animationing
for it = numel(t):-1:1
  for ip = 4:-1:1
    IePAD{ip}(:,:,it) = squeeze(Ie_alts{5-ip}([1:end,end],it,:))./...
                          (dE(1:size(IeZTE,3))'*BeamW([1:end,end]))';
  end
end

output = IePAD;
