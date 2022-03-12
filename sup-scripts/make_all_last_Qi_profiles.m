%% Volume and Column intensity-plots
% This script produces plots of altitude-time-variation of volume
% emission-rates and time-variation of normalized
% column-emission(excitation) intensity plots.

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

%% Optional setting of legend location:
if ~exist('leg_location')
  leg_location = 'northwestoutside';
end
% For possible options see documentation for legend.

%% Hard-coded output-directory
FigDir = sprintf('Figures-%s',datestr(now,'yyyymmdd'))
mkdir(fullfile(results_dir,FigDir))

ph = [];
legstr = {};
curr_idx = 1;
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
      try
        if ~isempty(strfind(dDir(i1).name,'DE'))
          C = 1/2;
        else
          C = 1;
        end
        load Qzt_all_L.mat  QOi QO2i QN2i h_atm t
        Qlast(:,curr_idx) = C*(QOi(:,end)+QN2i(:,end)+QO2i(:,end));
        ph(end+1) = plot(C*(QOi(:,end)+QN2i(:,end)+QO2i(:,end)),h_atm/1e3);
        hold on
        xlabel('ionization-rate (/m^3/s)','fontsize',15)
        ylabel('height (km)','fontsize',15)
        title('Ionization-rate','fontsize',15)
        ax = axis;
        %legstr{end+1} = dDir(i1).name;
        legstr{curr_idx} = dDir(i1).name;
        curr_idx = curr_idx + 1;
      catch
        fprintf('Everythings not right in directory: %s\n',pwd)
      end
    end
  end
end

legend(ph,...
       legstr,...
       'location',leg_location)


orient portrait
% Figures saved as both .png and .eps-files
print('-depsc2','-painters',fullfile(results_dir,FigDir,[dDir(i1).name,'-',RunDirs{i2},'-It-01.eps']))
print('-dpng','-painters',fullfile(results_dir,FigDir,[dDir(i1).name,'-',RunDirs{i2},'-It-01.png']))

%   Copyright � 2019 Bjorn Gustavsson, <bjorn.gustavsson@irf.se>
%   This is free software, licensed under GNU GPL version 2 or later
