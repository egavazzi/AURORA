
Science_grade_not_demo = 1;

if Science_grade_not_demo
  n_dirs = 1 + 180*4
  stepE = 1
else % Then demo with coarser energy-grid
  n_dirs = 1 + 180*4;
  stepE = 30;
end

%% Time intervall
t = linspace(0,0.025,31);


if ~exist('setup_completed','var') ||  setup_completed ~= 1
   % was this: setup4etrptd10streamsD
   print_figs = 0;
   plot_figs = 1;
   setup4etrptd10streamsB % setup4etrptdms
end

Isteady = 0e12;
I_flick = 0e12;

%% Runs with magnetic mirroring
par_list4DD = {350e3,   1,   3,       6, 1e3, 'I',  'I0_wBM-B6-350-1keV-1.0';
               350e3,   1,   3,       7, 1e3, 'I',  'I0_wBM-B7-350-1keV-1.0';
               350e3,   1,   3,       5, 1e3, 'I',  'I0_wBM-B5-350-1keV-1.0'};

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
i2do1 = 1:2:size(par_list4DD,1);
i2do2 = 2:2:size(par_list4DD,1);
i2do3 = size(par_list4DD,1):-1:1;
%%
for i_pars = i2do3, %Run_section(1)
  
  disp('===========================================================================')
  disp(par_list4DD(i_pars,:))
  disp('---------------------------------------------------------------------------')
  z_source = par_list4DD{i_pars,1}; % altitude of initial pulse of electron flux
  D_sc     = par_list4DD{i_pars,2}; % Diffusion coefficient scaling
  n_loop   = par_list4DD{i_pars,3};% number of time-steps to loop
  Beam_nr  = par_list4DD{i_pars,4};% 
  E0       = par_list4DD{i_pars,5}; % 
  mod_type = par_list4DD{i_pars,6}; % 'S' square-wave, 'H' sine-wave 'P' pulse, 'I' initial flux
  savedir  = par_list4DD{i_pars,7}; % '0Hz_H8keV';
  curr_par = par_list4DD(i_pars,:);
  
  %% Source altitude above ionosphere
  h_over_zmax = z_source - h_atm(end);
  ops4IeMstd.CscD_e = D_sc
  %% 0 Hertz flickering 8 keV
  % I0 = 0*squeeze(I00);
  [SUCCESS,MESSAGE,MESSAGEID] = mkdir(savedir);
  save(fullfile(savedir,'neutral_atm.mat'),'Te','h_atm','nN2','nO','nO2','ne')
  save(fullfile(savedir,'curr_par.mat'),'curr_par')
  
  
  %% The FAB-electron transport calculations
  %  First we have to specify the electron spectra. Peticolas and Lummerzheim
  %  used an electron spectra that had a constant number-flux per eV from the
  %  a peak energy at 3 keV all the way down to "the lowest energy" (which
  %  is not explicitly stated?). Here we chose to use 3 keV and 100 eV as
  %  energy-limits for the field-aligned burst. Further since we use
  %  multi-stream calculations we set the flux to be within 3 degrees from
  %  field-aligned.
  %  Peak-Energy of the FAB
  [~,iEu] = min(abs(E0-E))
  %  Bottom energy of the FAB
  [~,iEl] = min(abs(100-E));
  iEmax = iEu;
  
  % zeros in all directions
  for inots = numel(c_o_mu):-1:1,
    Ie_oftG{inots} = @(t,Ei) zeros(size(t));
  end
  % Except Primary electron-spectra FABs
  % flux only in a couple of field-aligned directions
  % Initial condition
  
  I0 = zeros(numel(h_atm)*(numel(mu_lims)-1),iEu);
  % and the electron-transport calculations
  [~,idxZ] = min(abs(z_source - h_atm))
  I0(numel(h_atm)*(Beam_nr-1)+idxZ,iEu) = 1e12;
  %% Loop over the time-intervall n_loop-times to split the problem
  %  into manageable chunks.
  for iN = 1:n_loop,
    fprintf('%d: %s\n',iN,datestr(now,'HH:MM:SS'))
    t_run = t(1:31)+(iN-1)*t(31);
    [Ie_ztE] = Ie_Mstream_tz_Bm_2_aurora_faster(h_atm,mag_ze,B,...
                                                E(1:iEmax),...
                                                mu_lims,mu_scatterings,...
                                                t_run,...
                                                I0,Ie_oftG,...
                                                zeros(numel(h_atm),numel(E(1:iEu))),...
                                                ne,Te,...
                                                ops4IeMstd,...
                                                nO,O_levels,XsO(:,1:iEmax),@O_e_2nd_dist,@phase_fcn_O,...
                                                nN2,N2_levels,XsN2(:,1:iEmax),@N2_e_2nd_dist,@phase_fcn_N2,...
                                                nO2,O2_levels,XsO2(:,1:iEmax),@O2_e_2nd_dist,@phase_fcn_O2);
    I0 = squeeze(Ie_ztE(:,end,:));
    savefile = fullfile(savedir,sprintf('IeFlickering-%02d.mat',iN));
    save(savefile,'Ie_ztE','E','t_run','mu_lims','h_atm','mu_scatterings','-v7.3')
  end
  
end


%% Runs WITHOUT magnetic mirroring
par_list4DD = {350e3,   1,   3,       6, 1e3, 'I',  'I0_nBM-B6-350-1keV-1.0';
               350e3,   1,   3,       7, 1e3, 'I',  'I0_nBM-B7-350-1keV-1.0';
               350e3,   1,   3,       5, 1e3, 'I',  'I0_nBM-B5-350-1keV-1.0'};

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
i2do1 = 1:2:size(par_list4DD,1);
i2do2 = 2:2:size(par_list4DD,1);
i2do3 = size(par_list4DD,1):-1:1;
%%
for i_pars = i2do3, %Run_section(1)
  
  disp('===========================================================================')
  disp(par_list4DD(i_pars,:))
  disp('---------------------------------------------------------------------------')
  z_source = par_list4DD{i_pars,1}; % altitude of initial pulse of electron flux
  D_sc     = par_list4DD{i_pars,2}; % Diffusion coefficient scaling
  n_loop   = par_list4DD{i_pars,3};% number of time-steps to loop
  Beam_nr  = par_list4DD{i_pars,4};% 
  E0       = par_list4DD{i_pars,5}; % 
  mod_type = par_list4DD{i_pars,6}; % 'S' square-wave, 'H' sine-wave 'P' pulse, 'I' initial flux
  savedir  = par_list4DD{i_pars,7}; % '0Hz_H8keV';
  curr_par = par_list4DD(i_pars,:);
  
  %% Source altitude above ionosphere
  h_over_zmax = z_source - h_atm(end);
  ops4IeMstd.CscD_e = D_sc
  %% 0 Hertz flickering 8 keV
  % I0 = 0*squeeze(I00);
  [SUCCESS,MESSAGE,MESSAGEID] = mkdir(savedir);
  save(fullfile(savedir,'neutral_atm.mat'),'Te','h_atm','nN2','nO','nO2','ne')
  save(fullfile(savedir,'curr_par.mat'),'curr_par')
  
  
  %% The FAB-electron transport calculations
  %  First we have to specify the electron spectra. Peticolas and Lummerzheim
  %  used an electron spectra that had a constant number-flux per eV from the
  %  a peak energy at 3 keV all the way down to "the lowest energy" (which
  %  is not explicitly stated?). Here we chose to use 3 keV and 100 eV as
  %  energy-limits for the field-aligned burst. Further since we use
  %  multi-stream calculations we set the flux to be within 3 degrees from
  %  field-aligned.
  %  Peak-Energy of the FAB
  [~,iEu] = min(abs(E0-E))
  %  Bottom energy of the FAB
  [~,iEl] = min(abs(100-E));
  iEmax = iEu;
  % zeros in all directions
  for inots = numel(c_o_mu):-1:1,
    Ie_oftG{inots} = @(t,Ei) zeros(size(t));
  end
  % Except Primary electron-spectra FABs
  % flux only in a couple of field-aligned directions
  % Initial condition
  
  I0 = zeros(numel(h_atm)*(numel(mu_lims)-1),iEu);
  % and the electron-transport calculations
  [~,idxZ] = min(abs(z_source - h_atm))
  I0(numel(h_atm)*(Beam_nr-1)+idxZ,iEu) = 1e12;
  %% Loop over the time-intervall n_loop-times to split the problem
  %  into manageable chunks.
  for iN = 1:n_loop,
    fprintf('%d: %s\n',iN,datestr(now,'HH:MM:SS'))
    t_run = t(1:31)+(iN-1)*t(31);
    [Ie_ztE] = Ie_Mstream_tz_Bm_2_aurora_faster(h_atm,mag_ze,ones(size(B)),...
                                                E(1:iEmax),...
                                                mu_lims,mu_scatterings,...
                                                t_run,...
                                                I0,Ie_oftG,...
                                                zeros(numel(h_atm),numel(E(1:iEu))),...
                                                ne,Te,...
                                                ops4IeMstd,...
                                                nO,O_levels,XsO(:,1:iEmax),@O_e_2nd_dist,@phase_fcn_O,...
                                                nN2,N2_levels,XsN2(:,1:iEmax),@N2_e_2nd_dist,@phase_fcn_N2,...
                                                nO2,O2_levels,XsO2(:,1:iEmax),@O2_e_2nd_dist,@phase_fcn_O2);
    I0 = squeeze(Ie_ztE(:,end,:));
    savefile = fullfile(savedir,sprintf('IeFlickering-%02d.mat',iN));
    save(savefile,'Ie_ztE','E','t_run','mu_lims','h_atm','mu_scatterings','-v7.3')
  end
  
end



%% Lower energies higer z0
%% Runs with magnetic mirroring
par_list4DD = {400e3,   1,   2,       6, 1e2, 'I',  'I0_wBM-B6-400-100eV-1.0';
               400e3,   1,   2,       7, 1e2, 'I',  'I0_wBM-B7-400-100eV-1.0';
               400e3,   1,   2,       8, 1e2, 'I',  'I0_wBM-B8-400-100eV-1.0'
               350e3,   1,   2,       6, 1e1, 'I',  'I0_wBM-B6-350-10eV-1.0';
               350e3,   1,   2,       7, 1e1, 'I',  'I0_wBM-B7-350-10eV-1.0';
               350e3,   1,   2,       8, 1e1, 'I',  'I0_wBM-B8-350-10eV-1.0'
              };
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
i2do1 = 1:2:size(par_list4DD,1);
i2do2 = 2:2:size(par_list4DD,1);
i2do3 = size(par_list4DD,1):-1:1;
%%
for i_pars = i2do3, %Run_section(1)
  
  disp('===========================================================================')
  disp(par_list4DD(i_pars,:))
  disp('---------------------------------------------------------------------------')
  z_source = par_list4DD{i_pars,1}; % altitude of initial pulse of electron flux
  D_sc     = par_list4DD{i_pars,2}; % Diffusion coefficient scaling
  n_loop   = par_list4DD{i_pars,3};% number of time-steps to loop
  Beam_nr  = par_list4DD{i_pars,4};% 
  E0       = par_list4DD{i_pars,5}; % 
  mod_type = par_list4DD{i_pars,6}; % 'S' square-wave, 'H' sine-wave 'P' pulse, 'I' initial flux
  savedir  = par_list4DD{i_pars,7}; % '0Hz_H8keV';
  curr_par = par_list4DD(i_pars,:);
  
  %% Source altitude above ionosphere
  h_over_zmax = z_source - h_atm(end);
  ops4IeMstd.CscD_e = D_sc
  %% 0 Hertz flickering 8 keV
  % I0 = 0*squeeze(I00);
  [SUCCESS,MESSAGE,MESSAGEID] = mkdir(savedir);
  save(fullfile(savedir,'neutral_atm.mat'),'Te','h_atm','nN2','nO','nO2','ne')
  save(fullfile(savedir,'curr_par.mat'),'curr_par')
  
  
  %% The FAB-electron transport calculations
  %  First we have to specify the electron spectra. Peticolas and Lummerzheim
  %  used an electron spectra that had a constant number-flux per eV from the
  %  a peak energy at 3 keV all the way down to "the lowest energy" (which
  %  is not explicitly stated?). Here we chose to use 3 keV and 100 eV as
  %  energy-limits for the field-aligned burst. Further since we use
  %  multi-stream calculations we set the flux to be within 3 degrees from
  %  field-aligned.
  %  Peak-Energy of the FAB
  [~,iEu] = min(abs(E0-E))
  %  Bottom energy of the FAB
  [~,iEl] = min(abs(100-E));
  iEmax = iEu;
  
  % zeros in all directions
  for inots = numel(c_o_mu):-1:1,
    Ie_oftG{inots} = @(t,Ei) zeros(size(t));
  end
  % Except Primary electron-spectra FABs
  % flux only in a couple of field-aligned directions
  % Initial condition
  
  I0 = zeros(numel(h_atm)*(numel(mu_lims)-1),iEu);
  % and the electron-transport calculations
  [~,idxZ] = min(abs(z_source - h_atm))
  I0(numel(h_atm)*(Beam_nr-1)+idxZ,iEu) = 1e12;
  %% Loop over the time-intervall n_loop-times to split the problem
  %  into manageable chunks.
  for iN = 1:n_loop,
    fprintf('%d: %s\n',iN,datestr(now,'HH:MM:SS'))
    t_run = t(1:31)+(iN-1)*t(31);
    [Ie_ztE] = Ie_Mstream_tz_Bm_2_aurora_faster(h_atm,mag_ze,B,...
                                                E(1:iEmax),...
                                                mu_lims,mu_scatterings,...
                                                t_run,...
                                                I0,Ie_oftG,...
                                                zeros(numel(h_atm),numel(E(1:iEu))),...
                                                ne,Te,...
                                                ops4IeMstd,...
                                                nO,O_levels,XsO(:,1:iEmax),@O_e_2nd_dist,@phase_fcn_O,...
                                                nN2,N2_levels,XsN2(:,1:iEmax),@N2_e_2nd_dist,@phase_fcn_N2,...
                                                nO2,O2_levels,XsO2(:,1:iEmax),@O2_e_2nd_dist,@phase_fcn_O2);
    I0 = squeeze(Ie_ztE(:,end,:));
    savefile = fullfile(savedir,sprintf('IeFlickering-%02d.mat',iN));
    save(savefile,'Ie_ztE','E','t_run','mu_lims','h_atm','mu_scatterings','-v7.3')
  end
  
end

%% Lower energies higer z0
%% Runs with magnetic mirroring
par_list4DD = {400e3,   1,   2,       6, 1e2, 'I',  'I0_nBM-B6-400-100eV-1.0';
               400e3,   1,   2,       7, 1e2, 'I',  'I0_nBM-B7-400-100eV-1.0';
               400e3,   1,   2,       8, 1e2, 'I',  'I0_nBM-B8-400-100eV-1.0'
               350e3,   1,   2,       6, 1e1, 'I',  'I0_nBM-B6-350-10eV-1.0';
               350e3,   1,   2,       7, 1e1, 'I',  'I0_nBM-B7-350-10eV-1.0';
               350e3,   1,   2,       8, 1e1, 'I',  'I0_nBM-B8-350-10eV-1.0'
              };
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
i2do1 = 1:2:size(par_list4DD,1);
i2do2 = 2:2:size(par_list4DD,1);
i2do3 = size(par_list4DD,1):-1:1;
%%
for i_pars = i2do3, %Run_section(1)
  
  disp('===========================================================================')
  disp(par_list4DD(i_pars,:))
  disp('---------------------------------------------------------------------------')
  z_source = par_list4DD{i_pars,1}; % altitude of initial pulse of electron flux
  D_sc     = par_list4DD{i_pars,2}; % Diffusion coefficient scaling
  n_loop   = par_list4DD{i_pars,3};% number of time-steps to loop
  Beam_nr  = par_list4DD{i_pars,4};% 
  E0       = par_list4DD{i_pars,5}; % 
  mod_type = par_list4DD{i_pars,6}; % 'S' square-wave, 'H' sine-wave 'P' pulse, 'I' initial flux
  savedir  = par_list4DD{i_pars,7}; % '0Hz_H8keV';
  curr_par = par_list4DD(i_pars,:);
  
  %% Source altitude above ionosphere
  h_over_zmax = z_source - h_atm(end);
  ops4IeMstd.CscD_e = D_sc
  %% 0 Hertz flickering 8 keV
  % I0 = 0*squeeze(I00);
  [SUCCESS,MESSAGE,MESSAGEID] = mkdir(savedir);
  save(fullfile(savedir,'neutral_atm.mat'),'Te','h_atm','nN2','nO','nO2','ne')
  save(fullfile(savedir,'curr_par.mat'),'curr_par')
  
  
  %% The FAB-electron transport calculations
  %  First we have to specify the electron spectra. Peticolas and Lummerzheim
  %  used an electron spectra that had a constant number-flux per eV from the
  %  a peak energy at 3 keV all the way down to "the lowest energy" (which
  %  is not explicitly stated?). Here we chose to use 3 keV and 100 eV as
  %  energy-limits for the field-aligned burst. Further since we use
  %  multi-stream calculations we set the flux to be within 3 degrees from
  %  field-aligned.
  %  Peak-Energy of the FAB
  [~,iEu] = min(abs(E0-E))
  %  Bottom energy of the FAB
  [~,iEl] = min(abs(100-E));
  iEmax = iEu;
  
  % zeros in all directions
  for inots = numel(c_o_mu):-1:1,
    Ie_oftG{inots} = @(t,Ei) zeros(size(t));
  end
  % Except Primary electron-spectra FABs
  % flux only in a couple of field-aligned directions
  % Initial condition
  
  I0 = zeros(numel(h_atm)*(numel(mu_lims)-1),iEu);
  % and the electron-transport calculations
  [~,idxZ] = min(abs(z_source - h_atm))
  I0(numel(h_atm)*(Beam_nr-1)+idxZ,iEu) = 1e12;
  %% Loop over the time-intervall n_loop-times to split the problem
  %  into manageable chunks.
  for iN = 1:n_loop,
    fprintf('%d: %s\n',iN,datestr(now,'HH:MM:SS'))
    t_run = t(1:31)+(iN-1)*t(31);
    [Ie_ztE] = Ie_Mstream_tz_Bm_2_aurora_faster(h_atm,mag_ze,ones(size(B)),...
                                                E(1:iEmax),...
                                                mu_lims,mu_scatterings,...
                                                t_run,...
                                                I0,Ie_oftG,...
                                                zeros(numel(h_atm),numel(E(1:iEu))),...
                                                ne,Te,...
                                                ops4IeMstd,...
                                                nO,O_levels,XsO(:,1:iEmax),@O_e_2nd_dist,@phase_fcn_O,...
                                                nN2,N2_levels,XsN2(:,1:iEmax),@N2_e_2nd_dist,@phase_fcn_N2,...
                                                nO2,O2_levels,XsO2(:,1:iEmax),@O2_e_2nd_dist,@phase_fcn_O2);
    I0 = squeeze(Ie_ztE(:,end,:));
    savefile = fullfile(savedir,sprintf('IeFlickering-%02d.mat',iN));
    save(savefile,'Ie_ztE','E','t_run','mu_lims','h_atm','mu_scatterings','-v7.3')
  end
  
end
