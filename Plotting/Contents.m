% PLOTTING
% Functions for animating and plotting the evolution of the
% time-varying energetic electron fluxes, Ie([h x mu],t,E)
% 
% time-energy-altitude plots of multi-stream electron fluxes
% 
% Animations of electron number fluxes
%  animate_IeztE_3DEzoft  - Ie(z,t,E) plotted as I(E,z) animated along t
%  animate_IeztE_3DtEofz  - Ie(z,t,E) plotted as I(t,E) animated along z
%  animate_IeztE_3DtzofE  - Ie(z,t,E) plotted as I(t,z) animated along E
%  animate_IeEztE_3DtlinEofz - Ie(z,t,E) plotted as I(t,E) animated along z
% Animations of electron energy fluxes
%  animate_IeEztE_3DEzoft - Ie(z,t,E)*E plotted as I(E,z) animated along t
%  animate_IeEztE_3DtEofz - Ie(z,t,E)*E plotted as I(t,E) animated along z
% Animations of pitch-angle distributions at selected altitudes
%  animate_IeztE_3DEthetaoft - pitch-angle electron energy flux animation
%  animate_IeztE_pitchangleflux - pitch-angle electron energy flux animation
% 
% Plots of electron number fluxes
%  plot_IeztE_3DEzoft     - Ie(z,t,E) plotted as I(E,z) animated along t
%  plot_IeztE_3DtEofz     - Ie(z,t,E) plotted as I(t,E) animated along z
%  plot_IeztE_3DtzofE     - Ie(z,t,E) plotted as I(t,z) animated along E
%  plot_IeztE_3DEthetaoft - pitch-angle plot of time-varying electron energy flux
%  plot_IezE_2DEtheta     - energy-pitch-angle electron flux plot
% Plots of electron energy fluxes
%  plot_IeEztE_3DEzoft - Ie(z,t,E)*E plotted as I(E,z) animated along t
%  plot_IeEztE_3DtEofz - Ie(z,t,E)*E plotted as I(t,E) animated along z
% 
% Plots of phase-space-density
%  plot_fezv_2Dvtheta - energy-pitch-angle phase-space-density plot
%
% Plots of excited states 
%  plot_exc1state      - UNTITLED2 Summary of this function goes here
%  plot_exc5states     - UNTITLED2 Summary of this function goes here
%  plot_excMstates     - UNTITLED2 Summary of this function goes here
%
% Figure-generating scripts:
%  Fig_PnL5n6        - time-energy plots of electron-spectra at 2 altitudes
%  Fig_PnL_Iezt_E    - time-altitude plots of electron fluxes at 4 energies
%  Fig_PnL_VEMCEM    - 4278 Å volume and column emission plot
%  threeDplot1stream - Script producing a 3-cuts plot

%   Copyright © 2019 Bjorn Gustavsson, <bjorn.gustavsson@irf.se>
%   This is free software, licensed under GNU GPL version 2 or later
