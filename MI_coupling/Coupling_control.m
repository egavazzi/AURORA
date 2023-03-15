%% This script controls the M-I coupling
% It calls AURORA and ketchup and converts the data between them.

%{
 -  Place yourself in the folder where AURORA is to be run, 
    i.e. /mnt/data/etienne/AURORA/MI_coupling/conversion_1.27e7_120s-2/
 -  Run part 1 (if not first run)
 -  Run part 2
 -  Run part 3
 -  Run part 4
 -  Run part 5 (if not first run)
 -  Run part 6
%}


%% To run from AURORA simulation folder
current_folder = pwd;
% Set path to the vlasov simulation state from which we are starting the coupling
path_to_vlasov_initial_file = "/mnt/data/etienne/ketchup/MI_coupling/l7/1.27e7/1.27e7_160to300s/outp/fzvzmu5600000.mat"; 
% Set path to the vlasov simulation to extracts from
path_to_vlasov_simulation = "/mnt/data/etienne/ketchup/MI_coupling/l7/1.27e7_MI-2";
path_to_next_vlasov_simulation = "/mnt/data/etienne/ketchup/MI_coupling/l7/1.27e7_MI-3";
path_to_vlasov_simulation_results = append(path_to_vlasov_simulation,'_results');

%% 1. Convert the .dat results from the ketchup simulation to .mat
if ~exist(path_to_vlasov_simulation_results)
  disp('Creating result folder ...')
  command = strcat("cp -r ",path_to_vlasov_simulation," ",path_to_vlasov_simulation_results);
  system(command);
  cd(path_to_vlasov_simulation_results)
  disp('... succesfully created!')
  disp('Starting the conversion of data from .dat to .mat ...')
  ketchup_b6conv
  disp('... done!')
  mkdir fzvzmu_boundary;
  movefile('outp/fzvzmuIB*.mat','fzvzmu_boundary/');
  mkdir fzvzmu;
  movefile('outp/fzvzmu*.mat','fzvzmu/');
  cd(current_folder)
end

%% 2. Convert the distr. funct. from KETCHUP to fluxes that can be used in AURORA
index_specie = 1;
HMR_VZ = 20; % HOW MUCH DO YOU WANT TO REFINE VZ
HMR_MU = 20; % HOW MUCH DO YOU WANT TO REFINE MU
E_max = 7231; % (eV) maximum energy to consider
theta_lims2do = [180:-10:90]; % pitch-angle limits for the streams (180° is straight downward)
first_run = 0; % if 1 it will extract only the initial state (constant incoming flux)

Ie_total = conversion_M_to_I(path_to_vlasov_initial_file,path_to_vlasov_simulation_results,...
                                index_specie,HMR_VZ,HMR_MU,E_max,theta_lims2do,first_run);
                                
disp('Saving the incoming fluxes ...')                                
save('Ie_incoming.mat','-v7.3','Ie_total') % saving in ketchup/
cd(current_folder)
save('Ie_incoming.mat','-v7.3','Ie_total') % saving in AURORA/
disp('... saved!')

%% 3. Calculate the ionospheric response
% Here we set the parameters
%                n_t,   Dir-name
par_list4G_10 = { 6,   'MIC-18streams-0.30s'};
Loop_MI_coupling

%% 4. Extract fluxes at the top of ionosphere
make_all_Ie_top_MI

%% 5. Prepare the folder for next vlasov simulation

if ~exist(path_to_next_vlasov_simulation)
  disp('Creating folder for next vlasov simulation ...')
  command = strcat("cp -r ",path_to_vlasov_simulation," ",path_to_next_vlasov_simulation);
  system(command);
  cd(path_to_next_vlasov_simulation)
  !rm dumps/*
  !find outp/datfiles/ -name "*.dat" -delete
  disp('... created!')
  cd(current_folder)
end

%% 6. Convert the fluxes from AURORA to distr. funct. that can be used in KETCHUP
index_specie = 1;
HMR_E = 5; % HOW MUCH DO YOU WANT TO REFINE E
HMR_MU = 20; % HOW MUCH DO YOU WANT TO REFINE MU
E_max = 7231; % (eV) maximum energy to consider
theta_lims2do = [180:-10:0]; % pitch-angle limits for the streams (180° is straight downward)

[f_out,f_out2] = conversion_I_to_M(path_to_vlasov_initial_file,path_to_vlasov_simulation_results, ...
                                      current_folder,index_specie,HMR_E,HMR_MU,E_max,theta_lims2do);

disp('Saving outgoing fluxes...')
% saving in AURORA/
fileID = fopen('f_response.bin','w');
fwrite(fileID',f_out,'double');   
fclose(fileID);
fileID = fopen('f_response2.bin','w');
fwrite(fileID',f_out2,'double');
fclose(fileID);
% saving in ketchup/
cd(path_to_next_vlasov_simulation)
fileID = fopen('f_response.bin','w');
fwrite(fileID',f_out,'double');   
fclose(fileID);
fileID = fopen('f_response2.bin','w');
fwrite(fileID',f_out2,'double');
fclose(fileID);
cd(current_folder)
disp('... saved!')