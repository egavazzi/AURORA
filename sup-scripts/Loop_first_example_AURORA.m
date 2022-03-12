%% Example-script for using AURORA the time-dependent electron-transport code
%  This shows an example usage of AURORA for calculating the
%  ionospheric fluxes for a short time-period with precipitation in
%  field-aligned bursts from 100 eV to 3 keV with 8 Hz
%  modulation-frequency a, from a source-altitude at 3000 km.

%% First set-up
% For science-grade we use a ~12 eV energy-grid (smaller grid-size
% at energies below 30 eV), and for the angular scattering we use
% 721 sub-angles to calculate the scattering matrices.
Science_grade_not_demo = 1;

if Science_grade_not_demo
  n_dirs = 1 + 180*4
  stepE = 1
else % Then demo with coarser energy-grid
  n_dirs = 1 + 180*4;
  stepE = 30
end

%% Time intervall
t = linspace(0,0.35,106);

%% Electron-fluxes
Isteady = 1e12;
I_flick = 1e12;

if ~exist('setup_completed','var') ||  setup_completed ~= 1
  setup4etrptd10streams
end

%% Specification of FAB-characteristics
% pars for HET:     z0,  f, n_t  BeamWeights E_max  mod  Dir-name
%                  (m)  Hz   #           %     eV  S/H/P
par_list4PnL = {3000e3,  8,  4,  [1 1]*B_W(1), 3e3, 'H','PnLh-08-3000-3keV-0-30';
               };

%% Initial conditions 
% The initial condition in all pitch-angle-streams is zero, that is
% we start with zero fluxes in all streams at the beginning of the 

I00 = zeros(numel(h_atm)*(numel(mu_lims)-1),numel(E));

% Energy grid for the low-energy tail
E_LET = E;
% Initially we use many runs without the LET, so simply set a
% zero-flux in all pitch-angle-streams:
for i1 = (numel(mu_lims)-1):-1:1,
  Ie_LET{i1} = 0*E_LET;
end



%% Time-dependent electron transport example
%  This script illustrates the use of time-dependent multi-stream
%  electron transport calculations by reproducing the calculations
%  from Peticolas and Lummerzheim [2000].
i2do = [1:size(par_list4PnL,1)];
% This is a suitable place to split the work into parts and run the
% calculations on multiple computers, if running the code in
% cell-mode just assign i2do2 or i2do3 to i2do1

%%
for i_pars = i2do, %Run_section(1)
  
  disp('===========================================================================')
  disp(par_list4PnL(i_pars,:))
  disp('---------------------------------------------------------------------------')
  curr_par = par_list4PnL(i_pars,:);
  z_source = par_list4PnL{i_pars,1} % 3000e3 4000e3; % Source altitude at 3000 km
  f_f      = par_list4PnL{i_pars,2} % 5, 10
  n_loop   = par_list4PnL{i_pars,3}; % 5, 7;
  curr_BW  = par_list4PnL{i_pars,4}; % B_W;
  Emax     = par_list4PnL{i_pars,5} % 
  mod_type = par_list4PnL{i_pars,6} % 'S' square-wave, 'H' sine-wave
  savedir  = par_list4PnL{i_pars,7} % '0Hz_H8keV';
  savedirstr = sprintf('PnL%s-%02d-%d-%dkeV-%d-%d',...
                       lower(mod_type),...
                       f_f,...
                       z_source/1e3,...
                       Emax/1e3,...
                       180-theta_lims2do(1),...
                       180-theta_lims2do(1+numel(curr_BW)))
  if isequal(savedir,savedirstr)
    disp('That was a close shave...')
  else
    disp('...but here you were trying to cut your own throat...')
    savedir = savedirstr;
  end
  
  %% Source altitude above ionosphere
  h_over_zmax = z_source - h_atm(end);
  
  %% 0 Hertz flickering 8 keV
  I0 = 0*squeeze(I00);
  [SUCCESS,MESSAGE,MESSAGEID] = mkdir(savedir);
  
  
  %% The FAB-electron transport calculations
  %  First we have to specify the electron spectra. Peticolas and Lummerzheim
  %  used an electron spectra that had a constant number-flux per eV from the
  %  a peak energy at 3 keV all the way down to "the lowest energy" (which
  %  is not explicitly stated?). Here we chose to use 3 keV and 100 eV as
  %  energy-limits for the field-aligned burst. Further since we use
  %  multi-stream calculations we set the flux to be within 3 degrees from
  %  field-aligned.
  %  Peak-Energy of the FAB
  [~,iEu] = min(abs(Emax-E));
  %  Bottom energy of the FAB
  [~,iEl] = min(abs(100-E));

  % zeros in all directions
  for inots = numel(c_o_mu):-1:1,
    Ie_oftG{inots} = @(t,Ei) zeros(size(t));
  end
  % Except Primary electron-spectra FABs
  % flux only in a couple of field-aligned directions
  i_Hots = 1:numel(curr_BW)
  for ihots = i_Hots,
    switch mod_type
     case 'S'
      Ie_oftG{ihots} = @(t,Ei) Ie_PnL_flickeringS(t,Ei,1e-2,...
                                                  E(iEl:iEu),dE(iEl:iEu),...
                                                  f_f*2*pi,...
                                                  h_over_zmax,...
                                                  abs(real(c_o_mu(ihots))),1)*curr_BW(ihots)/sum(curr_BW);
     otherwise
      Ie_oftG{ihots} = @(t,Ei) Ie_PnL_flickeringH(t,Ei,1e-2,...
                                                  E(iEl:iEu),dE(iEl:iEu),...
                                                  f_f*2*pi,...
                                                  h_over_zmax,...
                                                  abs(real(c_o_mu(ihots))),1)*curr_BW(ihots)/sum(curr_BW);
    end
  end
  % Initial condition
  I0 = zeros(numel(h_atm)*(numel(mu_lims)-1),1);
  % and the electron-transport calculations
  
  [SUCCESS,MESSAGE,MESSAGEID] = mkdir(savedir);
  save(fullfile(savedir,'neutral_atm.mat'),'Te','h_atm','nN2','nO','nO2','ne')
  save(fullfile(savedir,'curr_par.mat'),'curr_par')
  %% Loop over the time-intervall n_loop-times to split the problem
  %  into manageable chunks.
  for iN = 1:n_loop,
    fprintf('%d: %s\n',iN,datestr(now,'HH:MM:SS'))
    t_run = t(1:31)+(iN-1)*t(31);
    [Ie_ztE,mu_scatterings] = Ie_Mstream_tz_2_aurora(h_atm,mag_ze,...
                                                     E(1:iEu),...
                                                     mu_lims,mu_scatterings,t_run,...
                                                     I0,Ie_oftG,...
                                                     zeros(numel(h_atm),numel(E(1:iEu))),...
                                                     ne,Te,...
                                                     [],...
                                                     nO,O_levels,XsO(:,1:iEu),@O_e_2nd_dist,@phase_fcn_O,...
                                                     nN2,N2_levels,XsN2(:,1:iEu),@N2_e_2nd_dist,@phase_fcn_N2,...
                                                     nO2,O2_levels,XsO2(:,1:iEu),@O2_e_2nd_dist,@phase_fcn_O2);
    for iB = 1:numel(Ie_oftG),
      for iE = iEu:-1:1,
        for it = numel(t_run):-1:1,
          Ie_source(iB,it,iE) = Ie_PnL_flickeringH(t_run(it),E(iE),1e-2,...
                                                  E(iEl:iEu),dE(iEl:iEu),...
                                                  f_f*2*pi,...
                                                  0.1,... % I think
                                                  abs(real(c_o_mu(ihots))),1)*curr_BW(ihots)/sum(curr_BW);
        end
      end
    end
    I0 = squeeze(Ie_ztE(:,end,:));
    savefile = fullfile(savedir,sprintf('IeFlickering-%02d.mat',iN));
    save(savefile,'Ie_ztE','E','t_run','mu_lims','h_atm','mu_scatterings','-v7.3')
  end
  
end

%   Copyright © 2019 Bjorn Gustavsson, <bjorn.gustavsson@irf.se>
%   This is free software, licensed under GNU GPL version 2 or later

