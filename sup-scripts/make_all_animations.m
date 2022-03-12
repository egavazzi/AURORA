%% Make all animations-script
% This script attempts to produce 4 animations per
% results-directory, 3 with plane-cuts through the z-t-E
% electron-flux array for each pitch-angle-stream and one with
% pitch-angle-spectra at four altitudes.

%% Intensity-limits. 
%  To make it possible or easy to have consistent intensity-limits
%  across multiple result-files there is currently no automatic
%  intensity-limits set in this script. This means that an array
%  for the colour-limits, cxmax has to be generated manually. 

%% Root result-directories
if ~exist('results_dir','var') || isempty(results_dir)
  disp('Please enter a "results_dir"')
  return
end

%% Run-directories:
if ~exist('RunDirs','var') || isempty(RunDirs)
  disp('Please enter one or several "RunDirs"')
  return
end
if ~exist('movies2make')
  movies2make = [1 1 1 1 1];
end

%% This sub-plot layout might need to be modified to adapt to layouts
% more suitable to other pitch-angle-stream configurations:
spp = [2*ones(10,1),5*ones(10,1),[1:5,(10:-1:6)]'];

%% These figure-sizes might need to be modified to adapt to other
%pitch-angle-stream configurations too.
figPwide = [692, 556, 997, 418];   % trying this one for 10-beam calculations
figPsquare = [1021, 421, 668, 553];% trying this one for 9-beam
                                   % calculations and for
                                   % pitch-angle animation
figPHigh = [1214, 326, 475, 648];
fig_sz = figPwide;


figure
whitebg([1 1 1])


%% Loop away!
for i2 = 1:numel(RunDirs)
  cd(results_dir)
  cd(RunDirs{i2})
  dDir = dir;
  
  % Variable for diagnostics
  iFaulties = 1;
  %# cxmax = 12*ones(numel(dDir),2);
  for iRF = 3:numel(dDir)
    cd(results_dir)
    cd(RunDirs{i2})
    if dDir(iRF).isdir
      try
        % loading the electron-transport results
        [t,h_atm,E,mu_lims,IeZTE,mu_scatterings] = Ie_ztE_loader({dDir(iRF).name});
        dE = diff(E);dE(end+1) = dE(end);
        BeamW = 2*pi*mu_scatterings{3};
        for i1 = numel(BeamW):-1:1,
          theta_str{i1} = sprintf('%3.1f - %3.1f',...
                                  180-180/pi*acos(mu_lims(i1)),...
                                  180-180/pi*acos(mu_lims(i1+1)));
        end
        
        try
          % Producing the first animation with subplots for energy
          % (x) - altitude (y) with time-variation 
          if movies2make(1)
            filename = fullfile(dDir(iRF).name,'IeztE_3DEzoft.avi');
            fprintf('Making animation: \n',filename)
            colormap(jet)
            set(gcf,'position',fig_sz)
            animate_IeztE_3DEzoft(t,h_atm,E(1:size(IeZTE,3)),...
                                  IeZTE,...
                                  dE(1:size(IeZTE,3)),BeamW,...
                                  [-5 0]+max(cxmax(min(end,iRF),:)),spp, theta_str,filename);
          end
        catch
          ERR_MSG{iFaulties} = lasterr;
          disp(['something wrong with doing: ',filename])
          Faulties{iFaulties} = filename;
          iFaulties = iFaulties + 1;
        end
        try
          % Producing the second animation with subplots for energy
          % (x) - altitude (y) with energy-variation 
          if movies2make(2)
            filename = fullfile(dDir(iRF).name,'IeztE_3DtzofE.avi');
            fprintf('Making animation: \n',filename)
            colormap(jet)
            set(gcf,'position',fig_sz)
            animate_IeztE_3DtzofE(t,h_atm,E(1:size(IeZTE,3)),...
                                  IeZTE,...
                                  dE(1:size(IeZTE,3)),BeamW,...
                                  [-5 0]+max(cxmax(min(end,iRF),:)),spp, theta_str,filename);
          end
        catch
          ERR_MSG{iFaulties} = lasterr;
          disp(['something wrong with doing: ',filename])
          Faulties{iFaulties} = filename;
          iFaulties = iFaulties + 1;
        end
        try
          % Producing the third animation with subplots for time
          % (x) - Energy (y) with altitude-variation 
          if movies2make(3)
            filename = fullfile(dDir(iRF).name,'IeztE_3DtEofz.avi');
            fprintf('Making animation: \n',filename)
            set(gcf,'position',fig_sz)
            colormap(jet)
            animate_IeztE_3DtEofz(t,h_atm,E(1:size(IeZTE,3)),...
                                  IeZTE,...
                                  dE(1:size(IeZTE,3)),BeamW,...
                                  [-5 0]+max(cxmax(min(end,iRF),:)),spp, theta_str,filename);
          end
        catch
          ERR_MSG{iFaulties} = lasterr;
          disp(['something wrong with doing: ',filename])
          Faulties{iFaulties} = filename;
          iFaulties = iFaulties + 1;
        end
        try
          % Producing the fourth animation with pitch-angle
          % distribution at four heights
          if movies2make(4) 
            filename = fullfile(dDir(iRF).name,'IeztE_pitchangledist.avi');
            fprintf('Making animation: \n',filename)
            [dh,i115] = min(abs(h_atm/1e3-225));
            [dh,i175] = min(abs(h_atm/1e3-325));
            [dh,i300] = min(abs(h_atm/1e3-450));
            [dh,i600] = min(abs(h_atm/1e3-600));
            iZ = [i115,i175,i300,i600];
            set(gcf,'position',figPsquare)
            animate_IeztE_pitchangleflux(t,h_atm/1e3,...
                                         E,dE,...
                                         IeZTE,...
                                         BeamW,mu_lims,...
                                         iZ,...
                                         filename,...
                                         [-3.5 0]+max(cxmax(min(end,iRF),:)));
          end
        catch
          ERR_MSG{iFaulties} = lasterr;
          disp(['something wrong with doing: ',filename])
          Faulties{iFaulties} = filename;
          iFaulties = iFaulties + 1;
        end
        try
          % Producing the animation with fluxes as frunction of pitch-angle
          % and height at highest energy
          if movies2make(5)
            filename = fullfile(dDir(iRF).name,'IeztE_mu_z_at_E.avi');
            set(gcf,'position',figPHigh)
            fprintf('Making animation: \n',filename)
            animate_IeztE_3DzmuatEoft(t,h_atm,E,...
                                      IeZTE,...
                                      dE,BeamW,...
                                      size(IeZTE,3),...
                                      [6.5 12.5]-4,...
                                      {'0','10','30','60','80','90','100','120','150','170','180'},...
                                      filename);
          end
        catch
          ERR_MSG{iFaulties} = lasterr;
          disp(['something wrong with doing: ',filename])
          Faulties{iFaulties} = filename;
          iFaulties = iFaulties + 1;
        end
      catch
        ERR_MSG{iFaulties} = lasterr;
        disp(['something wrong with processing directory: ',dDir(iRF).name])
        Faulties{iFaulties} = dDir(iRF).name;
        iFaulties = iFaulties + 1;
      end
    end
  end
end


%   Copyright � 2019 Bjorn Gustavsson, <bjorn.gustavsson@irf.se>
%   This is free software, licensed under GNU GPL version 2 or later
