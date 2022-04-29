%% Time-varying electron fluxes for given incoming flux
% This script calculates the ionospheric electron-fluxes for a given
% input file containing fluxes coming down.

%% Scicence-grade run.
% Set this parameter to zero for a faster run to check that the
% script runs properly, with such poor energy-resolution that the
% results are physically useless.
Science_grade_not_demo = 1;

if Science_grade_not_demo
  n_dirs = 1 + 180*4
  stepE = 1
else % Then demo with coarser energy-grid
  n_dirs = 1 + 180*4;
  stepE = 30;
end

%% Time intervall
% The time-intervall per integration-period is deliberaltely
% short. For each set of parametres the integration of the
% electron-transport equations is repeated 6 times, for each step
% the time-period is extended by a factor of 2, continuing from the
% last time-step. This takes the solution to 3.15 s after the start.
t = 0:0.005:0.05;

%% Some arbitrary peak number of electron-flux
%Isteady = 1e12;
%I_flick = 1e12;

if ~exist('setup_completed','var') ||  setup_completed ~= 1
  % then we run the setup where neutral atmosphere, cross-sections
  % secondary-electron-spectra, phase-functions, scattering
  % matrices etc are calculated.
  setup4etrptd36streams % setup4etrptdms
end

%% Here we set the parameters
%                n_t,   Dir-name
%                       char
par_list4G_10 = {  3,   'MIC-36streams-0.35s--1'};


%% Let's go!
I00 = zeros(numel(h_atm)*(numel(mu_lims)-1),numel(E));

% Energy grid for the low-energy tail
E_LET = E;
% Initially we use many runs without the LET, so simply set a
% zero-flux in all pitch-angle-streams:
for i1 = (numel(mu_lims)-1):-1:1,
  Ie_LET{i1} = 0*E_LET;
end

%  This section calculates the rise-characteristics of
%  electron-fluxes in the ionosphere for mono-energetic
%  Precipitation, with energy and pitch-angle distribution as
%  parameterized in PAR_LIST4G_FA

  disp('===========================================================================')
  disp(AURORA)
  disp('---------------------------------------------------------------------------')
  n_loop   = par_list4G_10{i_pars,1}; % 3 loops --> 0.35s
  savedir  = par_list4G_10{i_pars,2}
  curr_par = par_list4G_10(i_pars,:)
  
  h_over_zmax = 0;
  iEmax = find(    cumsum(exp(-(E-E0dE(1)).^2/E0dE(2)^2)) < ...
               0.9999*sum(exp(-(E-E0dE(1)).^2/E0dE(2)^2)),1,'last')

  I0 = 0*squeeze(I00(:,1:iEmax));
  [SUCCESS,MESSAGE,MESSAGEID] = mkdir(savedir);
  save(fullfile(savedir,'neutral_atm.mat'),'Te','h_atm','nN2','nO','nO2','ne')
  save(fullfile(savedir,'curr_par.mat'),'curr_par')
  t_run = t;

  for i1 = 1:n_loop,
    fprintf('%d: %s\n',i1,datestr(now,'HH:MM:SS'))
    p_e_q = zeros(length(h_atm),length(E(1:iEmax)));
    
    [Ie_ztE] = Ie_Mstream_tz_2_aurora_MI(AURORA_root_directory,h_atm,mag_ze,E(1:iEmax),mu_lims,mu_scatterings,...
                                      t_run,...
                                      I0,p_e_q,ne,Te,...
                                      [],...
                                      nO,O_levels,XsO(:,1:iEmax),@O_e_2nd_dist,@phase_fcn_O,...
                                      nN2,N2_levels,XsN2(:,1:iEmax),@N2_e_2nd_dist,@phase_fcn_N2,...
                                      nO2,O2_levels,XsO2(:,1:iEmax),@O2_e_2nd_dist,@phase_fcn_O2);
    t_run = t_run(end) + 2*(t_run-t_run(1));
    I0 = squeeze(Ie_ztE(:,end,:));
    savefile = fullfile(savedir,sprintf('IeFlickering-%02d.mat',i1));
    save(savefile,'Ie_ztE','E','t_run','mu_lims','h_atm','Ie_oftGB','I0','mu_scatterings','-v7.3')
  end
  
end % Run_section(1)