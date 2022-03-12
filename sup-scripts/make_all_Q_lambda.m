%% Volume-emission-calculations
% This script calculates the volume-emission(excitation)-rates
% for the time-dependent electron-transport fluxes.
% This script batch-processes all results-directories in the
% RunDirs-cell-array with results-directories. The catual
% calculations are done with the Ie_ztE2Q_zt functions.

%% Root result-directories
if ~exist('results_dir','var') || isempty(results_dir)
  results_dir = '/media/bgu001/5f5e8978-a828-4fd4-aabf-2032a3fb895b/Data/Bjorns-panic-repository/Etrp-tdms';
end

%% Run-directories:
if ~exist('RunDirs','var') || isempty(RunDirs)
  RunDirs = {'Run_PnL20190828';
             'Run_AWA20190920';
             'Run_PnL20190916';
             'Run_PnL20191018'};
end

%% Loop away!
for i2 = 1:numel(RunDirs)
  cd(results_dir)
  cd(RunDirs{i2})
  dDir = dir;
  for i1 = 3:numel(dDir)
    cd(results_dir)
    cd(RunDirs{i2})
    if dDir(i1).isdir
      cd(dDir(i1).name)
      Ie_matfiles = dir('IeFlickering-*.mat');
      Qzt_file = dir('Qzt_all_L.mat');
      atm_file = dir('neutral_atm.mat');
      if ~isempty(atm_file)
        load(atm_file.name,'nN2','nO2','nO','h_atm')
      end
      its_done = 0;
      if ~isempty(Qzt_file)
        if all(Qzt_file.datenum - [Ie_matfiles.datenum] > 0)
          its_done = 1;
        end
      end
      if its_done
        CD = pwd;
        fprintf('Directory: %s already processed?\n',CD)
      else
        ThatWentOK = Ie_ztE2Q_zt(nN2,nO2,nO,h_atm);
        if ThatWentOK
          CD = pwd;
          fprintf(':::Processed OK: %s\n',CD)
        else
          CD = pwd;
          fprintf('Everythings not right in directory: %s\n',CD)
        end
      end
    end
  end
end


%   Copyright © 2019 Bjorn Gustavsson, <bjorn.gustavsson@irf.se>
%   This is free software, licensed under GNU GPL version 2 or later
