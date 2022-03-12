%  Utilities functions
%  Path ->  /home/bgu001/matlab/AURORA/utilities
%  Assorted aeronomy-specific functions and code-snippets
% 
% Functions for handling the time-dependent multi-stream electron fluxes
%   Ie_ztE2Q_zt      - volume emission/excitation rates, time-depentendly
%   Ie_zE2Qz         - volume emission/excitation rates, steady-state
%   Ie_pad2fe        - electron flux to phase-space density conversion
%   Ie_ztE_loader    - loads time-varying electron-flux data from AURORA
%   Ie_zE_loader     - loads steady-state electron-flux data from AURORA
%   exc_tz_of_Ie_ztE - volume excitation rate altitude-time variation
%   exc_z_of_Ie_zE   - steady-state volume excitation rate
%   Ie2currents      - calculate current from multistream electron fluxes
%   q2colem          - Integrate volume excitation-rates to column emission
%
% Aeronomy-related utility functions
%   E_levelsOnN2       - wavelength and enegy of 6 states in O and N2
%   e_meanfreepath     - electron mean-free-path
%   e_scattering_depth - electron scattering depth
%   get_all_xs         - electron-neutral impact cross sections collector-function
%   v_of_E             - convert electron energy (eV) to velocity (m/s)
%   E_SCATTERING_RESULT_FINDER - scattering-matrix-loading/calculating
%
% Signal-processing-function
%   dt_avgstd             - average time-shift in unsynchronised time-series
% 
% Plotting-functions
%   animate_auroral_rays  - Animate Auroral Rays
%   I_lambda_trajectories - stylish plot(I_1(t),I_2(t)) Useless function


%   Copyright © 2019 Bjorn Gustavsson, <bjorn.gustavsson@irf.se>
%   This is free software, licensed under GNU GPL version 2 or later
