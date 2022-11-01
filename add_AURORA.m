%% Set-up of matlab-path for the AURORA toolbox


%% Finding the path to the toolbox
[S1,S2,S3] = fileparts(which('add_AURORA'));

%% append the directory and sub-directories the path
% at the end! If some unfortunate name-clash occurs we do not want
% to break matlab's internal working, so accept that we have to fix
% our names to mend the name-clashes...
addpath(S1,'-end') % /home/bgu001/matlab/AURORA/,'-end')
addpath(fullfile(S1,'Data_electron'),'-end')
addpath(fullfile(S1,'Data_neutrals'),'-end')
addpath(fullfile(S1,'Data_photon'),'-end')
addpath(fullfile(S1,'Data_photon/N2'),'-end')
addpath(fullfile(S1,'Data_photon/NO'),'-end')
addpath(fullfile(S1,'Data_photon/O2'),'-end')
addpath(fullfile(S1,'Data_photon/O'),'-end')
addpath(fullfile(S1,'Demo-data'),'-end')    
addpath(fullfile(S1,'e_transport_Ms'),'-end')
addpath(fullfile(S1,'e_spectras'),'-end')   
addpath(fullfile(S1,'Plotting'),'-end')
addpath(fullfile(S1,'tools'),'-end') 
addpath(fullfile(S1,'tools/cm_and_cb_utilities'),'-end')
addpath(fullfile(S1,'tools/Factorize'),'-end')          
addpath(fullfile(S1,'tools/Geometry'),'-end') 
addpath(fullfile(S1,'utilities'),'-end') 
addpath(fullfile(S1,'tools/IGRF'),'-end') 
addpath(fullfile(S1,'sup-scripts'),'-end') 
addpath(fullfile(S1,'MI_coupling'),'-end')
addpath(fullfile(S1,'MI_coupling/src'),'-end')

AURORA_root_directory = S1;
%% Not to clutter workspace
clear S1 S2 S3
