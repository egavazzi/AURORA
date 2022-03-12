%% AERONOMY - Tools for auroral electron transport 
% aeronomy is a toolbox for solving time-dependent multi-stream electron
% transport equations for energetic electrons (ca 2 eV - 30 keV). The
% scripts this directory illustrates the use of the TDMS-e-transport code
% by calculating the electron-fluxes for several electron-spectra varying
% in time with 5 and 10 Hz.
%
% Example-scripts
%
%   Etrp_flickering_PnL                 - Time-dependent electron transport example
%   Etrp_4albedo8keV_9stream            - Steady-state albedo fluxes calculated with Time-dependent electron transport
%   Etrp_4LET8keV_9stream               - Steady-state Low-Energy-Tail for 8 keV Mono-Energetic Precipitation
%   Etrp_Flickering_8keVLET_9stream     - Steady-state 8 keV monoenergetic plus Low-Energy-Tail precipitation
%   Etrp_Flickering_HET_9stream         - Electron-fluxes for High-Energy Flickering precipitation
%   Etrp_Flickering_HETpLET8keV_9stream - Electron-fluxes for High-Energy Flickering precipitation
%   Etrp_Flickering_PnLpLET8keV_9stream - Electron fluxes for Low-Energy Field-aligned Bursts at 5 and 10 Hz
%   Etrp_Flickering_PnLpNothing_9stream - Calculation of steady-state Low-Energy-Tail 
%   Etrp_4LET2to8keV_9stream            - Steady-state Low-Energy-Tail for 2-8 keV mono-energetic precipitation
%   Etrp_Multi_streams                  - Time-dependent multi-stream electron transport example
%   setup4etrptdms                      - Set-up script for Examples in time-dependent electron transport
% 
% Industrial-strength run-it-all
%  
%   run_all_flickering                  - Run all the scripts above in
%   order
% 
% Result analysis scripts
%  
%   calc_all_flickering_emissions - Flickering excitation and ionization profiles
%   Flickering_tmp                      - 
%   pitch_angle_check                   - 
%   present_flickering_PnL              - Calculate excitation and ionization rates:
% 
% Aeronomy-functions
% 
%   dEds_ee          - suprathermal electron energy loss function for e-e-collisions 
%   dEdt_ee          - suprathermal electron energy loss rate for e-e-collisions 
%   get_all_xs       - collect all electron-neutral impact cross sections
%   v_of_E           - Transform electron energy (eV) to velocity (m/s)
%   exc_tz_of_Ie_ztE - calculate volume excitation or ionization rates
% 
% Path and documentation files
%   add_aeronomy                        - 
%   m2html4_aeronomy_2htmldoc           - Script run to generate cross-referenced html-documentation of code.
