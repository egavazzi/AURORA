function output = animate_IeztE_3DEthetaoft(t,h_atm,E,Ie_ztE,dE,BeamSA,mu_lims,cx_lims,spp, theta_strs,movieOut)
% animate_IeEztE_3DtzofE - animation of time-varying electron energy flux
% stepping in time. 
% animate_IeztE_3DEzoft animates the Energy and altitude-variation of
% electron energy-fluxes time-step by time-step from the start to the
% finish. The electron spectra are plotted in log-scale
% with pcolor in a number of subplots.
% 
% Calling:
%  clims = animate_IeztE_3DtEofz(t,h_atm,E,Ie_ztE,dE,BeamSA,cx_lims,spp, theta_strs)
% Input:
%  t          - time-scale (s), double array [1 x n_t]
%  h_atm      - altitudes (km), double array [n_z x 1]
%  E          - energies (eV), double array [1 x n_E]
%  Ie_ztE     - electron number-flux, double array [n_z x n_beams,n_t,n_E]
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
%  movieOut   - Optional argument, if MovieName in matlab-supported
%               video-format then animate_IeEztE_3DtEofz will attempt
%               to write the displayed sequence to a video-file
%               with that filename using the VideoWriter
%               functionality. See VideoWriter for details. If
%               numerical value that is converted to TRUE, a matlab
%               movie will be returned.
% Output:
%  clims - 3-D array with color-scales for all sub-plots plotted.
%  M     - Matlab-movie with the displayed sequence.
% 
% Example: See Etrp_example_7stream.m

%   Copyright © 2018 Björn Gustavsson, <bjorn.gustavsson@uit.no>
%   This is free software, licensed under GNU GPL version 2 or later

nZ = numel(h_atm);
doMovie = 0;

if nargin > 9 && ~isempty(movieOut)
  if isstr(movieOut)
    try
      vidObj = VideoWriter(movieOut);
      open(vidObj);
      doMovie = 1;
    catch
      doMovie = 0;
      disp(['Failed to initialize VideoWriter Object for movie: ',movieOut])
    end
  elseif movieOut
    doMovie = 2;
  end
end
iZ = [ 72   227   377];

for it = 1:numel(t),
  subplot(3,1,1)
  pcolor(v_of_E(E),180-180/pi*acos(mu_lims),log10(squeeze(Ie_ztE([iZ(3):numel(h_atm):end,end],it,:))./(BeamSA([1:end,end])'*dE))),shading flat
  set(gca,...
      'xtick',v_of_E([1 10 100 300 1000 3000 1e4]),...
      'xticklabel',{'1' '10' '100' '300','1000','3000', '10,000'},...
      'ytick',[0 30 60 90 120 150 180],...
      'tickdir','out',...
      'box','off')
  if isempty(cx_lims)
    caxis([-4 0]+max(caxis))
  else
    caxis(cx_lims)
  end
  title(sprintf('electron flux at %3.1f km, t = %3.3f s',h_atm(iZ(3))/1e3,t(it)),'fontsize',15)
  subplot(3,1,2)
  pcolor(v_of_E(E),180-180/pi*acos(mu_lims),log10(squeeze(Ie_ztE([iZ(2):numel(h_atm):end,end],it,:))./(BeamSA([1:end,end])'*dE))),shading flat
  set(gca,...
      'xtick',v_of_E([1 10 100 300 1000 3000 1e4]),...
      'xticklabel',{'1' '10' '100' '300','1000','3000', '10,000'},...
      'ytick',[0 30 60 90 120 150 180],...
      'tickdir','out',...
      'box','off')
  if isempty(cx_lims)
    caxis([-4 0]+max(caxis))
  else
    caxis(cx_lims)
  end
  ylabel('pitch-angle','fontsize',15)
  title(sprintf('electron flux at %3.1f km, t = %3.3f s',h_atm(iZ(2))/1e3,t(it)),'fontsize',15)
  colorbar_labeled('/eV/m^2/s/ster','log')
  subplot(3,1,3)
  pcolor(v_of_E(E),180-180/pi*acos(mu_lims),log10(squeeze(Ie_ztE([iZ(1):numel(h_atm):end,end],it,:))./(BeamSA([1:end,end])'*dE))),shading flat
  set(gca,...
      'xtick',v_of_E([1 10 100 300 1000 3000 1e4]),...
      'xticklabel',{'1' '10' '100' '300','1000','3000', '10,000'},...
      'ytick',[0 30 60 90 120 150 180],...
      'tickdir','out',...
      'box','off')
  if isempty(cx_lims)
    caxis([-4 0]+max(caxis))
  else
    caxis(cx_lims)
  end
  title(sprintf('electron flux at %3.1f km, t = %3.3f s',h_atm(iZ(1))/1e3,t(it)),'fontsize',15)
  xlabel('Energy (eV)','fontsize',15)
  drawnow
  if doMovie == 1 || doMovie == 2
    if doMovie == 1
      writeVideo(vidObj,getframe(gcf));
    elseif doMovie == 2
      M(it) = getframe(gcf);
    end
  end
end

if doMovie == 2
  output = M;
elseif doMovie == 1
  close(vidObj);
else
  output = [];
end
