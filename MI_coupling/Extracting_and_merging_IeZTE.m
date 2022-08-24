% This script can be used to extract electron fluxes from the ketchup data,
% and from the AURORA data. It then creates an animation of the electron
% fluxes as a heatmap over altitude and energy

current_folder = pwd;
% Set path to the vlasov simulation state from which we are starting the coupling
path_to_vlasov_initial_file = "/mnt/data/etienne/ketchup/MI_coupling/l7/1.27e7/1.27e7_160to300s/outp/fzvzmu5600000.mat"; 
% Set path to the vlasov simulation to extracts from
path_to_vlasov_simulation = "/mnt/data/etienne/ketchup/MI_coupling/l7/dt_2e-3/1.27e7_MI-1";
path_to_vlasov_simulation_results = append(path_to_vlasov_simulation,'_results');

%% Getting electron fluxes from ionosphere
% Extract the data
[t,h_atm_iono,E,mu_lims,IeZTE_iono,mu_scatterings] = Ie_ztE_loader({'./MIC-18streams-0.30s/'});

%% Getting electron fluxes from magnetosphere
tic
% Choose until which altitude electron fluxes should be extracted.
zmax_to_extract = 1000; % in (km)

HMR_VZ = 5;   % HOW MUCH DO YOU WANT TO REFINE VZ
HMR_MU = 5;   % HOW MUCH DO YOU WANT TO REFINE MU

% Extract the data
[IeZTE_DL,z_DL] = make_Ie_from_ketchup(path_to_vlasov_initial_file,path_to_vlasov_simulation_results, ...
                                        index_specie,zmax_to_extract,h_atm_iono,HMR_VZ,HMR_MU, ...
                                        size(IeZTE_iono,3),mu_lims);
                                      
savefile = strcat(['IeZTE_DL_',num2str(zmax_to_extract),'km.mat']);
cd(path_to_vlasov_simulation_results);
save(savefile,'IeZTE_DL','z_DL','-v7.3');
cd(current_folder)
save(savefile,'IeZTE_DL','z_DL','-v7.3');
toc

%% Merging both Ie into one
for i_mu = 1:18
  IeZTE_iono_plot(307*(i_mu-1) + (1:307),:,:) = IeZTE_iono(310*(i_mu-1) + (1:307),:,:);
end
h_atm_iono_plot = h_atm_iono(1:307);
% IeZTE_iono_plot = IeZTE_iono;
% h_atm_iono_plot = h_atm_iono;

IeZTE = zeros(size(IeZTE_DL,1)+ size(IeZTE_iono_plot,1),size(IeZTE_iono_plot,2),size(IeZTE_iono_plot,3));
z = [h_atm_iono_plot/1e3;z_DL.'];
% initialisation
% time_2_update_DL = 0; % in (s)
% k = 1;
for n_t = 1:size(IeZTE,2)
%   if abs(t(n_t)-time_2_update_DL) < eps
%       % the time grid for our electron fluxes in the magnetosphere is
%       % coarser, so this 'if' block tests if it is time for an update of our
%       % magnetospheric electron fluxes
%       IeZTE_DL_temp = IeZTE_DL(:,k,:);                        % extract the Ie
%       time_2_update_DL = time_2_update_DL + 0.0025;             % and set the time of the next update
%       k = k + 1;
%   end
  for i_mu = 1:18
    IeZTE(numel(z)*(i_mu-1) + (1:numel(h_atm_iono_plot)),n_t,:) = ...
        IeZTE_iono_plot(numel(h_atm_iono_plot)*(i_mu-1) + (1:numel(h_atm_iono_plot)),n_t,:);
    IeZTE(numel(z)*(i_mu-1) + numel(h_atm_iono_plot) + (1:numel(z_DL)),n_t,:) = ...
        IeZTE_DL(numel(z_DL)*(i_mu-1) + (1:numel(z_DL)),n_t,:);
  end
end

%% Make a fancy video
z_2_plot = z;
spp = [2*ones(8,1),4*ones(8,1),[1:4,(8:-1:5)]'];
% spp = [1*ones(1,1),1*ones(1,1),[1]'];
results_dir = AURORA_root_directory;
RunDirs = {'MI_coupling/a_new_hope_dt_2e-3/conversion_1.27e7_120s-1'};


figPwide = [692, 556, 997, 418];   % trying this one for 10-beam calculations
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
        IeZTE_2_plot = zeros(numel(z_2_plot)*(numel(theta_lims_2_plot)-1),size(IeZTE,2),size(IeZTE,3));
        try
          for i1 = 1:numel(mu_lims_plot)-1
            % find the index of the streams to sum
            index1 = find(abs(mu_lims - mu_lims_plot(i1)) < 0.001);
            index2 = find(abs(mu_lims - mu_lims_plot(i1+1)) < 0.001);
            % and sum them
            for i2 = index1:(index2-1)
              hindex = (1:numel(z_2_plot))+(i1-1)*numel(z_2_plot);
              hindex_2_sum = (1:numel(z_2_plot))+(i2-1)*numel(z_2_plot);
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
          filename = fullfile(dDir(iRF).name,'IeztE_3DEzoft_1000km.avi');
          fprintf('Making animation: \n',filename)
          colormap(jet)
          set(gcf,'position',fig_sz)
          set(gcf,'WindowState','maximized');
          animate_IeztE_3DEzoft(t,z_2_plot,E(1:size(IeZTE,3)),...
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
        disp('meh')
%         Faulties{iFaulties} = filename;
        iFaulties = iFaulties + 1;
      end
    end
  end
end


