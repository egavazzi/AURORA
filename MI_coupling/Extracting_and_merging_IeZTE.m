% This script can be used to extract electron fluxes from the ketchup data,
% and from the AURORA data. It then creates an animation of the electron
% fluxes as a heatmap over altitude and energy

current_folder = pwd;

% Set path to the Vlasov simulation state from which we are starting the coupling
path_to_vlasov_initial_file = "/mnt/data/etienne/ketchup/MI_coupling/l7/1.27e7/Nz_3000_150to300s_results/f_vlasov_raw/fzvzmu6000000.mat";

% Set path to the Vlasov simulation we are interested in
path_to_vlasov_simulation = "/mnt/data/etienne/ketchup/MI_coupling/l7/1.27e7_MI/0.3s_total-5s";
path_to_vlasov_simulation_results = append(path_to_vlasov_simulation,'_results'); % This is only for old simulations. New simulation have results in same folder.

% Set path to the Aurora simulation we are interested in
path_to_aurora_simulation = {'/mnt/data/etienne/AURORA/MI_coupling/TIME2PLAY/conversion_1.27e7-1/MIC-18streams-0.30s/'};

%% Getting electron fluxes in the ionosphere
% Load the data
[t,h_atm,E,mu_lims,IeZTE_iono,mu_scatterings] = Ie_ztE_loader(path_to_aurora_simulation);

%% Getting electron fluxes in the magnetosphere
tic
% Choose until which altitude electron fluxes should be extracted.
zmax_to_extract = 10000; % (km)
HMR_VZ = 5;   % HOW MUCH DO YOU WANT TO REFINE VZ
HMR_MU = 5;   % HOW MUCH DO YOU WANT TO REFINE MU



filename = strcat(['IeZTE_DL_',num2str(zmax_to_extract),'km.mat']);
fullfilename = fullfile(path_to_vlasov_simulation_results, filename); % This is only for old simulations. New simulation have results in same folder.

loaded = 0;
% try to load the data if they already exist
if exist(fullfile(path_to_vlasov_simulation_results, filename))       % This is only for old simulations. New simulation have results in same folder.
  load(fullfilename, 'HMR_VZ_file', 'HMR_MU_file')                    % This is only for NEW simulations.
  if (HMR_VZ_file == HMR_VZ) && (HMR_MU_file == HMR_MU_file)
    load(fullfilename)
    loaded = 1;
  end
end
% if not, extract them and save them
if loaded == 0
  [IeZTE_DL,z_DL] = make_Ie_from_ketchup(path_to_vlasov_initial_file, path_to_vlasov_simulation_results, ...
                                          1, zmax_to_extract, h_atm_iono, HMR_VZ, HMR_MU, ...
                                          size(IeZTE_iono,3), mu_lims);% This is only for old simulations. New simulation have results in same folder.
  cd(path_to_vlasov_simulation_results);                              % This is only for old simulations. New simulation have results in same folder.
  save(filename,'IeZTE_DL','z_DL', 'HMR_VZ_file', 'HMR_MU_file', '-v7.3');
  cd(path_to_aurora_simulation)
  save(filename,'IeZTE_DL','z_DL', 'HMR_VZ_file', 'HMR_MU_file', '-v7.3');
end
toc
%% Merging both Ie into one
sampling_t_iono = 10;
sampling_t_DL = 5;
sampling_z_DL = 5;



clear IeZTE_iono_plot
clear IeZTE_DL_plot
% resampling Ie from the ionosphere
h_atm_plot = h_atm(end-2); % getting rid of the two points overlapping with the ionosphere
t_plot = t(1:sampling_t_iono:end);
for i_mu = 1:length(mu_lims)-1
  IeZTE_iono_plot(numel(h_atm_plot)*(i_mu-1) + (1:numel(h_atm_plot)), :, :) = ...
                                            IeZTE_iono(numel(h_atm)*(i_mu-1) + (1:numel(h_atm_plot)), 1:sampling_t_iono:end, :);
end
% resampling Ie from the magnetosphere
z_DL_plot = z_DL(1:sampling_z_DL:end);
for i_mu = 1:length(mu_lims)-1
  IeZTE_DL_plot(numel(z_DL_plot)*(i_mu-1) + (1:numel(z_DL_plot)), :, :) = ...
                                            IeZTE_DL(numel(z_DL)*(i_mu-1) + (1:5:numel(z_DL)), 1:sampling_t_DL:end, :);
end
% merging everything
IeZTE = zeros(size(IeZTE_DL_plot, 1) + size(IeZTE_iono_plot, 1), size(IeZTE_iono_plot, 2), size(IeZTE_iono_plot, 3));
z = [h_atm_plot/1e3; z_DL_plot.'];
for n_t = 1:size(IeZTE,2)
  for i_mu = 1:length(mu_lims)-1
    IeZTE(numel(z)*(i_mu-1) + (1:numel(h_atm_plot)), n_t, :) = ...
                IeZTE_iono_plot(numel(h_atm_plot)*(i_mu-1) + (1:numel(h_atm_plot)), n_t, :);
    IeZTE(numel(z)*(i_mu-1) + numel(h_atm_plot) + (1:numel(z_DL_plot)), n_t, :) = ...
                IeZTE_DL_plot(numel(z_DL_plot)*(i_mu-1) + (1:numel(z_DL_plot)), n_t, :);
  end
end

%% Only Ie in Magnetosphere
% clear IeZTE_iono_plot
% clear IeZTE_DL_plot
% for i_mu = 1:18
%   IeZTE_iono_plot(307*(i_mu-1) + (1:307),:,:) = IeZTE_iono(310*(i_mu-1) + (1:307),1:10:end,:);
% end
% h_atm_plot = h_atm_iono(1:307);
% z_DL_plot = z_DL(1:5:end);
% for i_mu = 1:18
%   IeZTE_DL_plot(numel(z_DL_plot)*(i_mu-1) + (1:numel(z_DL_plot)),:,:) = ...
%                                             IeZTE_DL(numel(z_DL)*(i_mu-1) + (1:5:numel(z_DL)),1:end,:);
% end
% 
% IeZTE = zeros(size(IeZTE_DL_plot,1),size(IeZTE_DL_plot,2),size(IeZTE_DL_plot,3));
% z = [z_DL_plot.'];
% 
% IeZTE = IeZTE_DL_plot;
% 
% t_plot = linspace(0, 5, size(IeZTE,2));
%% Make a fancy video
results_dir = AURORA_root_directory;
RunDirs = {'MI_coupling/TIME2PLAY/conversion_1.27e7-1'};
name_videofile = strcat('IeztE_3DEzoft_', num2str(zmax_to_extract),'_km.avi');


spp = [2*ones(8,1),4*ones(8,1),[1:4,(8:-1:5)]'];
figPwide = [710, 284, 1200, 610];   % trying this one for 10-beam calculations
fig_sz = figPwide;
figure
whitebg([1 1 1])

for ii = 1:numel(RunDirs)
  cd(results_dir)
  cd(RunDirs{ii})
  dDir = dir;
  
  % Variable for diagnostics
  iFaulties = 1;
  %# cxmax = 12*ones(numel(dDir),2);
  for iRF = 3:numel(dDir)
    cd(results_dir)
    cd(RunDirs{ii})
    if dDir(iRF).isdir
      try
        dE = diff(E);dE(end+1) = dE(end);
        BeamW = 2*pi*mu_scatterings{3}; % because sum of mu_scattering should be = 4pi
        theta_lims_2_plot = [180 160 130 110 90 70 50 20 0];
        % make a string for the plot titles, based on the given array
        % theta_lims_2_plot
        mu_lims_plot = cosd(theta_lims_2_plot);
        for i1 = numel(mu_lims_plot)-1:-1:1
          if mu_lims_plot(i1) < 0
            theta_str{i1} = sprintf('%3.1f - %3.1f DOWN',...
                          180-180/pi*acos(mu_lims_plot(i1)),...
                          180-180/pi*acos(mu_lims_plot(i1+1)));
          else
            theta_str{i1} = sprintf('%3.1f - %3.1f UP',...
                          180/pi*acos(mu_lims_plot(i1+1)),...
                          180/pi*acos(mu_lims_plot(i1)));
          end
        end
        % Sum the fluxes from the streams that are contained in between the
        % theta_lims_2_plot values.
        % Exemple: if we want to plot the total flux of electrons with
        % pitch-angles between 20° and 40°, it will sum the fluxes from the
        % streams between 20°-30° and 30°-40°, if these streams exist.
        IeZTE_2_plot = zeros(numel(z)*(numel(theta_lims_2_plot)-1),size(IeZTE,2),size(IeZTE,3));
        try
          for i1 = 1:numel(mu_lims_plot)-1
            % find the index of the streams to sum
            index1 = find(abs(mu_lims - mu_lims_plot(i1)) < 0.001);
            index2 = find(abs(mu_lims - mu_lims_plot(i1+1)) < 0.001);
            % and sum them
            for i2 = index1:(index2-1)
              hindex = (1:numel(z))+(i1-1)*numel(z);
              hindex_2_sum = (1:numel(z))+(i2-1)*numel(z);
              IeZTE_2_plot(hindex,:,:) =  IeZTE_2_plot(hindex,:,:) + IeZTE(hindex_2_sum,:,:)./BeamW(i2);
            end
          end
        catch
          disp(['Error : the pitch-angle to plot ', num2str(theta_lims_2_plot(i1)),...
             '° does not match any of the stream limits used in the simulation'])
          break
        end

        
        
        try
        % Producing the first animation with subplots for energy
        % (x) - altitude (y) with time-variation 
          filename = fullfile(dDir(iRF).name,name_videofile);
          fprintf('Making animation: \n',filename)
          colormap(jet)
          set(gcf,'position',fig_sz)
%           set(gcf,'WindowState','maximized');
          animate_IeztE_3DEzoft(t_plot,z,E(1:size(IeZTE,3)),...
                                IeZTE_2_plot,...
                                dE(1:size(IeZTE,3)),BeamW,...
                                [-7 2],spp, theta_str,filename);
        catch
          ERR_MSG{iFaulties} = lasterr;
          disp(['something wrong with doing: ',filename])
          Faulties{iFaulties} = filename;
          iFaulties = iFaulties + 1;
        end
      catch
        ERR_MSG{iFaulties} = lasterr;
        iFaulties = iFaulties + 1;
      end
    end
  end
end