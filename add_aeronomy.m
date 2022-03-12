%% Set-up of matlab-path for the aeronomy toolbox

[S1,S2,S3] = fileparts(which('add_aeronomy'));

addpath(S1,'-end') % /home/bgu001/matlab/aeronomy/,'-end')
addpath(fullfile(S1,'Data_electron'),'-end')
addpath(fullfile(S1,'Data_neutrals'),'-end')
addpath(fullfile(S1,'Demo-data'),'-end')    
addpath(fullfile(S1,'e_transport_Ms'),'-end')
addpath(fullfile(S1,'e_transport_2s'),'-end')
addpath(fullfile(S1,'e_spectras'),'-end')   
% addpath(fullfile(S1,'Flickering_examples'),'-end')
% addpath(fullfile(S1,'Steady_state_examples'),'-end')
addpath(fullfile(S1,'Plotting'),'-end')
addpath(fullfile(S1,'tools'),'-end') 
addpath(fullfile(S1,'tools/cm_and_cb_utilities'),'-end')
addpath(fullfile(S1,'tools/Factorize'),'-end')          
addpath(fullfile(S1,'tools/Geometry'),'-end') 
addpath(fullfile(S1,'utilities'),'-end') 
addpath(fullfile(S1,'tests'),'-end') 
addpath(fullfile(S1,'sup-scripts'),'-end') 

% $$$ addpath /home/bgu001/matlab/aeronomy/ -end
% $$$ addpath /home/bgu001/matlab/aeronomy/Data_electron/ -end
% $$$ addpath /home/bgu001/matlab/aeronomy/Data_neutrals/ -end
% $$$ addpath /home/bgu001/matlab/aeronomy/Demo-data/ -end    
% $$$ addpath /home/bgu001/matlab/aeronomy/e_transport_Ms/ -end
% $$$ addpath /home/bgu001/matlab/aeronomy/e_transport_2s/ -end
% $$$ addpath /home/bgu001/matlab/aeronomy/e_spectras/ -end    
% $$$ addpath /home/bgu001/matlab/aeronomy/Flickering_examples/ -end
% $$$ addpath /home/bgu001/matlab/aeronomy/Steady_state_examples/ -end
% $$$ addpath /home/bgu001/matlab/aeronomy/Plotting/ -end
% $$$ addpath /home/bgu001/matlab/aeronomy/tools/ -end   
% $$$ addpath /home/bgu001/matlab/aeronomy/tools/cm_and_cb_utilities/ -end
% $$$ addpath /home/bgu001/matlab/aeronomy/tools/Factorize/ -end          
% $$$ addpath /home/bgu001/matlab/aeronomy/tools/Geometry/ -end 
% $$$ addpath /home/bgu001/matlab/aeronomy/tests -end 
% $$$ addpath /home/bgu001/matlab/aeronomy/sup-scripts -end 
