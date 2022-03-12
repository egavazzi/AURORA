function output = plot_IeEztE_3DtEofz(t,h_atm,E,Ie_ztE,dE,BeamSA,cx_lims,spp, theta_strs,idx_z)
% plot_IeEztE_3DtzofE - animation of time-varying electron energy flux
% stepping in time. 
% plot_IeztE_3DEzoft animates the Energy and altitude-variation of
% electron energy-fluxes time-step by time-step from the start to the
% finish. The electron spectra are plotted in log-scale
% with pcolor in a number of subplots.
% 
% Calling:
%  clims = plot_IeEztE_3DtEofz(t,h_atm,E,Ie_ztE,dE,BeamSA,cx_lims,spp, theta_strs,idx_z)
%  plot_IeEztE_3DtEofz(t,h_atm,E,Ie_ztE,dE,BeamSA,cx_lims,spp, theta_strs,idx_z)
% Input:
%  t,         - time-scale (s), double array [1 x n_t]
%  h_atm      - altitudes (km), double array [n_z x 1]
%  E          - energies (eV), double array [1 x n_E]
%  Ie_ztE     - electron number-flux, double array [n_z,n_t,n_E]
%  dE         - energy bin sizes (eV), double array [1 x n_E]
%  BeamSA     - solid angle sizes of the beams, [n_beams x 1]
%  cx_lims    - limit for colour-scale double array [1 x 2],
%               optional argument, if left empty the fluxes of each
%               beam will be scaled thusly: caxis([-4 0]+max(caxis))
%  spp        - sub-plot-position, integer array [n_beams x 3] with
%               sub-plot layout and subplot position row-by-row,
%               the function is designed for two rows of sub-plots
%               with the idea that downward fluxes should be
%               plotted in the top row and the upward in the bottom
%               row.
%  theta_strs - titles for the subplots, cell-array with
%               pitch-angle-boundaries, has to be at least as long
%               as the number of sub-plots plotted.
%  idx_z      - index in altitude-dimension to plot, scalar integer in
%               range of <1, 2,... ,numel(h_atm) >
% Output:
%  clims - 3-D array with color-scales for all sub-plots plotted.
% 
% Example: See Etrp_example_7stream.m

%   Copyright © 2018 Björn Gustavsson, <bjorn.gustavsson@uit.no>
%   This is free software, licensed under GNU GPL version 2 or later

nZ = numel(h_atm);
clims = [];

for i_z = idx_z,
  clf
  for i1 = 1:size(spp,1),     
    
    % log10(max(0,real(squeeze(Ie_ztE(i_z+(i1-1)*numel(h_atm),:,:))./(ones(size(h_atm))*dE)/BeamSA(i1))))),
    sph(i1) = subplot(spp(i1,1),spp(i1,2),spp(i1,3));
    pcolor(t,...
           E,...
           log10(max(0,real(squeeze(Ie_ztE(i_z+(i1-1)*numel(h_atm),:,:))'./repmat(dE(:)./(E(:)+dE(:)/2),size(t))/BeamSA(i1))))),
    shading flat
    if isempty(cx_lims)
      clims(i_z,i1,:) = caxis;
      caxis([-4 0]+max(caxis))
    else
      caxis(cx_lims)
    end
    set(gca,'tickdir','out','yscale','log','ytick',[1,10,1000,10000,1e5])
    if mod(i1,spp(1,2))==0 || i1 == size(spp,1)
      colorbar_labeled('log10(#eV/m2/s/eV/ster)')
    elseif isempty(cx_lims)
      colorbar_labeled('')      
    end
    if mod(i1,spp(1,2))==1 %i1 == 1 || i1 == 4
      ylabel('energy (eV)')
    else
      set(gca,'yticklabel','')        
    end
    if i1 > spp(1,2)
      xlabel('time (s)')
    else
      set(gca,'xticklabel','')
    end
    if i1 == 2
      title(sprintf('%4.3f (km)',h_atm(i_z)/1e3))
    else
      title(sprintf('pitch-angle ~ %s',theta_strs{i1}))
    end
  end
  drawnow
end

linkaxes(sph)
