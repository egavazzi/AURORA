%% Time-varying electron fluxes for given incoming flux
% This script calculates the ionospheric electron-fluxes for a given
% input file containing fluxes coming down.

%% Time intervall
% The time-intervall per integration-period is deliberaltely
% short. For each set of parametres the integration of the
% electron-transport equations is repeated 6 times, for each step
% the time-period is extended by a factor of 2, continuing from the
% last time-step. This takes the solution to 3.15 s after the start.
t = 0:1e-3:0.05;
% t = 0:2e-3:0.10;
% 
% if ~exist('setup_completed','var') ||  setup_completed ~= 1
  % then we run the setup where neutral atmosphere, cross-sections
  % secondary-electron-spectra, phase-functions, scattering
  % matrices etc are calculated.
  setup4etrptd18streams
% end

%% Here we set the parameters
%                n_t,   Dir-name
%                       char
% par_list4G_10 = {  3,   'MIC-18streams-0.35s-1'};

%% Let's go!
I00 = zeros(numel(h_atm)*(numel(mu_lims)-1),numel(E));

%  This section calculates the rise-characteristics of
%  electron-fluxes in the ionosphere for mono-energetic
%  Precipitation, with energy and pitch-angle distribution as
%  parameterized in PAR_LIST4G_FA

disp('======================')
disp('======= AURORA =======')
disp('----------------------')
n_loop   = par_list4G_10{1}; % 3 loops --> 0.35s
savedir  = par_list4G_10{2}
curr_par = par_list4G_10(:);

h_over_zmax = 0;
E_max = 7231; % (eV) maximum energy to consider

[~,iE_max] = min(abs(E-E_max)); 
I0 = 0*squeeze(I00(:,1:iE_max));
[SUCCESS,MESSAGE,MESSAGEID] = mkdir(savedir);
save(fullfile(savedir,'neutral_atm.mat'),'Te','h_atm','nN2','nO','nO2','ne')
save(fullfile(savedir,'curr_par.mat'),'curr_par')
t_run = t;

for i1 = 1:n_loop,
  fprintf('%d: %s\n',i1,datestr(now,'HH:MM:SS'))
  p_e_q = zeros(length(h_atm),length(E(1:iE_max)));

  [Ie_ztE] = Ie_Mstream_tz_2_aurora_MI(AURORA_root_directory,i1,h_atm,mag_ze,E(1:iE_max),mu_lims,mu_scatterings,...
                                    t_run,...
                                    I0,p_e_q,ne,Te,...
                                    [],...
                                    nO,O_levels,XsO(:,1:iE_max),@O_e_2nd_dist,@phase_fcn_O,...
                                    nN2,N2_levels,XsN2(:,1:iE_max),@N2_e_2nd_dist,@phase_fcn_N2,...
                                    nO2,O2_levels,XsO2(:,1:iE_max),@O2_e_2nd_dist,@phase_fcn_O2);
                                  
  I0 = squeeze(Ie_ztE(:,end,:));
  savefile = fullfile(savedir,sprintf('IeFlickering-%02d.mat',i1));
  save(savefile,'Ie_ztE','E','t_run','mu_lims','h_atm','I0','mu_scatterings','-v7.3')
  t_run = t_run(end) + (t_run-t_run(1));
end