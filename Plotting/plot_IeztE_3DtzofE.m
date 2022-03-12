function H = plot_IeztE_3DtzofE(t,h_atm,E,Ie_ztE,dE,BeamSA,cx_lims,spp, theta_strs,idx_E)
% plot_IeztE_3DtzofE - animation of time-varying electron flux
% energy-by-energy. 
% plot_IeztE_3DtzofE animates the altitude and time-variation of
% electron number-fluxes energy-by-energy from the highest energy
% down to the lowest. The electron spectra are plotted in log-scale
% with pcolor in a number of subplots.
% 
% Calling:
%  clims = plot_IeztE_3DtzofE(t,h_atm,E,Ie_ztE,dE,BeamSA,cx_lims,spp, theta_strs,idx_E)
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
%  idx_E      - index in altitude-dimension to plot, scalar integer in
%               range of <1, 2,... ,numel(E) >
% Output:
%   H - array of pcolor-handles
% Example: See Etrp_example_7stream.m

%   Copyright © 2018 Björn Gustavsson, <bjorn.gustavsson@uit.no>
%   This is free software, licensed under GNU GPL version 2 or later

nZ = numel(h_atm);
clims = [];

for iE = idx_E
  
  clf
  for i1 = 1:size(spp,1),
    
    sph(i1) = subplot(spp(i1,1),spp(i1,2),spp(i1,3));
    h(i1) = pcolor(t,...
                   h_atm,...
                   log10(max(0,real(Ie_ztE((1:nZ)+(i1-1)*nZ,:,iE)/dE(iE)/BeamSA(i1)))));
    shading flat,
    % caxis([-5 0]+max(caxis)),colorbar_labeled('')
    if isempty(cx_lims)
      clims(iE,i1,:) = caxis;
      caxis([-4 0]+max(caxis))
    else
      caxis(cx_lims)
    end
    set(gca,'tickdir','out')
    if mod(spp(i1,3),spp(1,2))==0 % mod(i1,spp(1,2))==0 % || i1 == size(spp,1)
      if all(BeamSA == 1)
        colorbar_labeled('log10(#/m2/s/eV)')
      else
        colorbar_labeled('log10(#/m2/s/eV/ster)')
      end
    elseif isempty(cx_lims)
      colorbar_labeled('','log')
    end
    if mod(spp(i1,3),spp(1,2))==1 % mod(i1,spp(1,2))==1 %i1 == 1 || i1 == 4
      ylabel('height (km)')
    else
      set(gca,'yticklabel','')
    end
    if spp(i1,3) > (spp(1,2)*(spp(1,1)-1)) % i1 > (spp(1,2)*(spp(1,1)-1))
      xlabel('time (s)')
    else
      set(gca,'xticklabel','')
    end
    if spp(i1,3) == ceil(spp(1,2)/2) % was i1 == 2
                              % if i1 == 2
      if E(iE) > 1e3
        title(sprintf('%2.2f (keV)',E(iE)/1e3))
      else
        title(sprintf('%2.2f (eV)',E(iE)))
      end
    else
      % title(sprintf('pitch-angle ~ %s',theta_strs{i1}))
      %title(sprintf('\theta_B ~ %s',theta_strs{i1}))
      title(['\theta_B ~ ',theta_strs{i1}])
    end
  end
  drawnow
end

if nargout
  H = h;
end

linkaxes(sph)