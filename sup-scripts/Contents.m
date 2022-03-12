%%  Batch-processing scripts
% Directory with scripts for running the time-dependent electron
% transport code with purpose in a reproducible manner. The
% directory contains example scripts for flickering aurora and
% other 
% short-time precipitation events. In addition a couple of scripts
% for handling analysis of the resulting data files to produce
% ionization and excitation profiles and generate plots and
% animations of the results.
%
% Scripts for running the electron transport for different events
% 
% Loop_first_example_AURORA      - flickering F-A STEB-precipitation single parameters
% Loop_Run_PnL_flickering10beams - flickering F-A STEB-precipitation varying parameters
%
% Scripts for runnning the post-processing of the electron
% transport output
% 
% make_all_Q_lambda        - multi-directory excitation and ionization-rates calculations
% make_all_I_lambda        - multi-directory column-emission calculations
% make_all_IQ_lambda_plots - multi-directory volume/column-emission plots
% make_all_Ie_top          - multi-directory primary-precipitation-spectra extraction
%
% Scripts for generating plots and animations
% 
% make_all_animations        - multi-directory electron-flux animation production
% combine_n_extract_I_dt_Mm5 - intensity-modulation and time-shift plots
% pars_of_ItMm_dirnames.txt  - parameters of flickering-precipitation,
%                              see combine_n_extract_I_dt_Mm5 for
%                              how it is used and necessary.

%   Copyright © 2019 Bjorn Gustavsson, <bjorn.gustavsson@irf.se>
%   This is free software, licensed under GNU GPL version 2 or later
