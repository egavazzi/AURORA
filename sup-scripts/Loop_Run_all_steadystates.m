
Science_grade_not_demo = 1;

if Science_grade_not_demo
  n_dirs = 1 + 180*4
  stepE = 1
else % Then demo with coarser energy-grid
  n_dirs = 1 + 180*4;
  stepE = 30;
end

%% Time intervall
t = linspace(0,0.35,106);

Isteady = 1e12;
I_flick = 1e12;

if ~exist('setup_completed','var') ||  setup_completed ~= 1
  setup4etrptd10streams % setup4etrptdms
end
% TODO: E_LET and Ie_LET needs to be fixed
% pars for HET:     z0,     f, n_t     BeamWeights    E0dE   Dir-name
%                  (m)     Hz      #           %          eV   char
par_list4G_fa = {3000e3, 0.0001,   1,      B_W(1),  [7e3 100], 'G_fa_ss-7keV';
                 3000e3, 0.0001,   1,      B_W(1),  [5e3 100], 'G_fa_ss-5keV';
                 3000e3, 0.0001,   1,      B_W(1),  [3e3 100], 'G_fa_ss-3keV'
                 3000e3, 0.0001,   1,      B_W(1:5),  [7e3 100], 'G_iso_ss-7keV';
                 3000e3, 0.0001,   1,      B_W(1:5),  [5e3 100], 'G_iso_ss-5keV';
                 3000e3, 0.0001,   1,      B_W(1:5),  [3e3 100], 'G_iso_ss-3keV'};
I00 = zeros(numel(h_atm)*(numel(mu_lims)-1),numel(E));

% Energy grid for the low-energy tail
E_LET = E;
% Initially we use many runs without the LET, so simply set a
% zero-flux in all pitch-angle-streams:
for i1 = (numel(mu_lims)-1):-1:1,
  Ie_LET{i1} = 0*E_LET;
end

%% Found a suitable Gaussian flux at higher energy
%E0dE = 8*[1,0.1]*1e3;
%EfdEf = [E0dE(1)+E0dE(2)*1,E0dE(2)*4.374];

%% High-Energy steady-state precipitation 3000 km, 0-15 deg, 8.8 keV
%
%  This section calculates the rise-characteristics of
%  electron-fluxes in the ionosphere for 8.8 keV mono-energetic
%  field-aligned precipitation (0-15 degrees pitch-angle).
% %% Etrp_Flickering_HET4_9stream.m:savedir = 'Onset_H8p8keV


for i_pars = 1:size(par_list4G_fa,1) %Run_section(1)
  disp('===========================================================================')
  disp(par_list4G_fa(i_pars,:))
  disp('---------------------------------------------------------------------------')
  z_source = par_list4G_fa{i_pars,1};  % 2390e3 3390e3; % Source altitude at 3000 km
  f_f      = par_list4G_fa{i_pars,2}; % 0, 5, 10
  n_loop   = par_list4G_fa{i_pars,3}; % 5, 7;
  curr_BW  = par_list4G_fa{i_pars,4}; % B_W;
  E0dE     = par_list4G_fa{i_pars,5}  % ;
  savedir  = par_list4G_fa{i_pars,6}  % 'Onset_H8keV';
  curr_par = par_list4G_fa(i_pars,:)
  
  h_over_zmax = z_source - h_atm(end);
  %curr_BW = curr_BW.*(1/3).^(0:(numel(curr_BW)-1))
  iEmax = find(    cumsum(exp(-(E-E0dE(1)).^2/E0dE(2)^2)) < ...
               0.9999*sum(exp(-(E-E0dE(1)).^2/E0dE(2)^2)),1,'last')
  %% 0 Hertz flickering 8 keV
  [SUCCESS,MESSAGE,MESSAGEID] = mkdir(savedir);
  save(fullfile(savedir,'neutral_atm.mat'),'Te','h_atm','nN2','nO','nO2','ne')
  save(fullfile(savedir,'curr_par.mat'),'curr_par')
  
  for i1 = 1:n_loop,
    
    fprintf('%d: %s\n',i1,datestr(now,'HH:MM:SS'))
    
    p_e_q = zeros(length(h_atm),length(E(1:iEmax)));
    % disp([E0dE(1,1),E0dE(1,2)])
    
    for iBeam = numel(c_o_mu):-1:1
      Ie_oftGB{iBeam} = @(E) zeros(size(E));
    end
    i_Hots = 1:numel(curr_BW);
    for ihots = i_Hots,
      Ie_oftGB{1} = @(Ei) 1e12*exp(-(Ei-E0dE(1,1)).^2/E0dE(1,2))*curr_BW(ihots)/sum(curr_BW);
    end
    
    [Ie_zE] = Ie_M_stream_4_aurora(h_atm,mag_ze,E(1:iEmax),mu_lims,mu_scatterings,Ie_oftGB,p_e_q,ne,Te,...
                                   nO,O_levels,XsO(:,1:iEmax),@O_e_2nd_dist,@phase_fcn_O,...
                                   nN2,N2_levels,XsN2(:,1:iEmax),@N2_e_2nd_dist,@phase_fcn_N2,...
                                   nO2,O2_levels,XsO2(:,1:iEmax),@O2_e_2nd_dist,@phase_fcn_O2);
    
    savefile = fullfile(savedir,sprintf('IeSteady_1-%02d.mat',i1));
    save(savefile,'Ie_zE','E','mu_lims','h_atm','Ie_oftGB','mu_scatterings','-v7.3')
  end
  
  
end % Run_section(1)




%% 12-stream mu-limits:
dtheta = [0 10 15 20 20 15 10 10 15 20 20 15 10];
theta_lims2do = 180-cumsum(dtheta);
mu_lims = cos(theta_lims2do*pi/180);
plot_everyone = 0;
%% Calculate the beam-to-beam scattering and weigthings
%  This can be done in the function Ie_Mstream_tz_2_aurora, as was done
%  above where the mu_scatterings argument was set to a place-holder 1, or
%  before as done here. This is necessary to get the proper weights for
%  calculating the fluxes in different beams to more for example an
%  isotropic flux
n_dirs = 721;
[Pmu2mup,theta2beamW,BeamW] = e_scattering_beamdistribution(mu_lims,n_dirs);
theta_lims = theta_lims2do;
mu_scatterings = {Pmu2mup,theta2beamW,BeamW};
c_o_mu = mu_avg(mu_lims);
save e_s_b_12streams.mat theta_lims mu_lims Pmu2mup theta2beamW BeamW c_o_mu

for i_pars = 1:size(par_list4G_fa,1) %Run_section(1)
  disp('===========================================================================')
  disp(par_list4G_fa(i_pars,:))
  disp('---------------------------------------------------------------------------')
  z_source = par_list4G_fa{i_pars,1};  % 2390e3 3390e3; % Source altitude at 3000 km
  f_f      = par_list4G_fa{i_pars,2}; % 0, 5, 10
  n_loop   = par_list4G_fa{i_pars,3}; % 5, 7;
  curr_BW  = par_list4G_fa{i_pars,4}; % B_W;
  E0dE     = par_list4G_fa{i_pars,5}  % ;
  savedir  = par_list4G_fa{i_pars,6}  % 'Onset_H8keV';
  curr_par = par_list4G_fa(i_pars,:)
  
  h_over_zmax = z_source - h_atm(end);
  %curr_BW = curr_BW.*(1/3).^(0:(numel(curr_BW)-1))
  iEmax = find(    cumsum(exp(-(E-E0dE(1)).^2/E0dE(2)^2)) < ...
               0.9999*sum(exp(-(E-E0dE(1)).^2/E0dE(2)^2)),1,'last')
  %% 0 Hertz flickering 8 keV
  [SUCCESS,MESSAGE,MESSAGEID] = mkdir(savedir);
  save(fullfile(savedir,'neutral_atm.mat'),'Te','h_atm','nN2','nO','nO2','ne')
  save(fullfile(savedir,'curr_par.mat'),'curr_par')
  
  for i1 = 1:n_loop,
    
    fprintf('%d: %s\n',i1,datestr(now,'HH:MM:SS'))
    
    p_e_q = zeros(length(h_atm),length(E(1:iEmax)));
    % disp([E0dE(1,1),E0dE(1,2)])
    
    for iBeam = numel(c_o_mu):-1:1
      Ie_oftGB{iBeam} = @(E) zeros(size(E));
    end
    i_Hots = 1:numel(curr_BW);
    for ihots = i_Hots,
      Ie_oftGB{1} = @(Ei) 1e12*exp(-(Ei-E0dE(1,1)).^2/E0dE(1,2))*curr_BW(ihots)/sum(curr_BW);
    end
    
    [Ie_zE] = Ie_M_stream_4_aurora(h_atm,mag_ze,E(1:iEmax),mu_lims,mu_scatterings,Ie_oftGB,p_e_q,ne,Te,...
                                   nO,O_levels,XsO(:,1:iEmax),@O_e_2nd_dist,@phase_fcn_O,...
                                   nN2,N2_levels,XsN2(:,1:iEmax),@N2_e_2nd_dist,@phase_fcn_N2,...
                                   nO2,O2_levels,XsO2(:,1:iEmax),@O2_e_2nd_dist,@phase_fcn_O2);
    
    savefile = fullfile(savedir,sprintf('IeSteady_2-%02d.mat',i1));
    save(savefile,'Ie_zE','E','mu_lims','h_atm','Ie_oftGB','mu_scatterings','-v7.3')
  end
  
  
end % Run_section(1)


%% 18-stream mu-limits:
dtheta = [0 3 5 7 15 30 15 7 5 3 3 5 7 15 30 15 7 5 3]
theta_lims2do = 180-cumsum(dtheta);
mu_lims = cos(theta_lims2do*pi/180);
plot_everyone = 0;
%% Calculate the beam-to-beam scattering and weigthings
%  This can be done in the function Ie_Mstream_tz_2_aurora, as was done
%  above where the mu_scatterings argument was set to a place-holder 1, or
%  before as done here. This is necessary to get the proper weights for
%  calculating the fluxes in different beams to more for example an
%  isotropic flux
n_dirs = 721;
[Pmu2mup,theta2beamW,BeamW] = e_scattering_beamdistribution(mu_lims,n_dirs);
theta_lims = theta_lims2do;
mu_scatterings = {Pmu2mup,theta2beamW,BeamW};
c_o_mu = mu_avg(mu_lims);
save e_s_b_18streams.mat theta_lims mu_lims Pmu2mup theta2beamW BeamW c_o_mu


for i_pars = 1:size(par_list4G_fa,1) %Run_section(1)
  disp('===========================================================================')
  disp(par_list4G_fa(i_pars,:))
  disp('---------------------------------------------------------------------------')
  z_source = par_list4G_fa{i_pars,1};  % 2390e3 3390e3; % Source altitude at 3000 km
  f_f      = par_list4G_fa{i_pars,2}; % 0, 5, 10
  n_loop   = par_list4G_fa{i_pars,3}; % 5, 7;
  curr_BW  = par_list4G_fa{i_pars,4}; % B_W;
  E0dE     = par_list4G_fa{i_pars,5}  % ;
  savedir  = par_list4G_fa{i_pars,6}  % 'Onset_H8keV';
  curr_par = par_list4G_fa(i_pars,:)
  
  h_over_zmax = z_source - h_atm(end);
  %curr_BW = curr_BW.*(1/3).^(0:(numel(curr_BW)-1))
  iEmax = find(    cumsum(exp(-(E-E0dE(1)).^2/E0dE(2)^2)) < ...
               0.9999*sum(exp(-(E-E0dE(1)).^2/E0dE(2)^2)),1,'last')
  %% 0 Hertz flickering 8 keV
  [SUCCESS,MESSAGE,MESSAGEID] = mkdir(savedir);
  save(fullfile(savedir,'neutral_atm.mat'),'Te','h_atm','nN2','nO','nO2','ne')
  save(fullfile(savedir,'curr_par.mat'),'curr_par')
  
  for i1 = 1:n_loop,
    
    fprintf('%d: %s\n',i1,datestr(now,'HH:MM:SS'))
    
    p_e_q = zeros(length(h_atm),length(E(1:iEmax)));
    % disp([E0dE(1,1),E0dE(1,2)])
    
    for iBeam = numel(c_o_mu):-1:1
      Ie_oftGB{iBeam} = @(E) zeros(size(E));
    end
    i_Hots = 1:numel(curr_BW);
    for ihots = i_Hots,
      Ie_oftGB{1} = @(Ei) 1e12*exp(-(Ei-E0dE(1,1)).^2/E0dE(1,2))*curr_BW(ihots)/sum(curr_BW);
    end
    
    [Ie_zE] = Ie_M_stream_4_aurora(h_atm,mag_ze,E(1:iEmax),mu_lims,mu_scatterings,Ie_oftGB,p_e_q,ne,Te,...
                                   nO,O_levels,XsO(:,1:iEmax),@O_e_2nd_dist,@phase_fcn_O,...
                                   nN2,N2_levels,XsN2(:,1:iEmax),@N2_e_2nd_dist,@phase_fcn_N2,...
                                   nO2,O2_levels,XsO2(:,1:iEmax),@O2_e_2nd_dist,@phase_fcn_O2);
    
    savefile = fullfile(savedir,sprintf('IeSteady_3-%02d.mat',i1));
    save(savefile,'Ie_zE','E','mu_lims','h_atm','Ie_oftGB','mu_scatterings','-v7.3')
  end
  
  
end % Run_section(1)



%% 20-stream mu-limits:
dtheta = [0 3 5 7 10 20 20 10 7 5 3 3 5 7 10 20 20 10 7 5 3 ]
theta_lims2do = 180-cumsum(dtheta);
mu_lims = cos(theta_lims2do*pi/180);
plot_everyone = 0;
%% Calculate the beam-to-beam scattering and weigthings
%  This can be done in the function Ie_Mstream_tz_2_aurora, as was done
%  above where the mu_scatterings argument was set to a place-holder 1, or
%  before as done here. This is necessary to get the proper weights for
%  calculating the fluxes in different beams to more for example an
%  isotropic flux
n_dirs = 721;
[Pmu2mup,theta2beamW,BeamW] = e_scattering_beamdistribution(mu_lims,n_dirs);
theta_lims = theta_lims2do;
mu_scatterings = {Pmu2mup,theta2beamW,BeamW};
c_o_mu = mu_avg(mu_lims);
save e_s_b_20streams.mat theta_lims mu_lims Pmu2mup theta2beamW BeamW c_o_mu

for i_pars = 1:size(par_list4G_fa,1) %Run_section(1)
  disp('===========================================================================')
  disp(par_list4G_fa(i_pars,:))
  disp('---------------------------------------------------------------------------')
  z_source = par_list4G_fa{i_pars,1};  % 2390e3 3390e3; % Source altitude at 3000 km
  f_f      = par_list4G_fa{i_pars,2}; % 0, 5, 10
  n_loop   = par_list4G_fa{i_pars,3}; % 5, 7;
  curr_BW  = par_list4G_fa{i_pars,4}; % B_W;
  E0dE     = par_list4G_fa{i_pars,5}  % ;
  savedir  = par_list4G_fa{i_pars,6}  % 'Onset_H8keV';
  curr_par = par_list4G_fa(i_pars,:)
  
  h_over_zmax = z_source - h_atm(end);
  %curr_BW = curr_BW.*(1/3).^(0:(numel(curr_BW)-1))
  iEmax = find(    cumsum(exp(-(E-E0dE(1)).^2/E0dE(2)^2)) < ...
               0.9999*sum(exp(-(E-E0dE(1)).^2/E0dE(2)^2)),1,'last')
  %% 0 Hertz flickering 8 keV
  [SUCCESS,MESSAGE,MESSAGEID] = mkdir(savedir);
  save(fullfile(savedir,'neutral_atm.mat'),'Te','h_atm','nN2','nO','nO2','ne')
  save(fullfile(savedir,'curr_par.mat'),'curr_par')
  
  for i1 = 1:n_loop,
    
    fprintf('%d: %s\n',i1,datestr(now,'HH:MM:SS'))
    
    p_e_q = zeros(length(h_atm),length(E(1:iEmax)));
    % disp([E0dE(1,1),E0dE(1,2)])
    
    for iBeam = numel(c_o_mu):-1:1
      Ie_oftGB{iBeam} = @(E) zeros(size(E));
    end
    i_Hots = 1:numel(curr_BW);
    for ihots = i_Hots,
      Ie_oftGB{1} = @(Ei) 1e12*exp(-(Ei-E0dE(1,1)).^2/E0dE(1,2))*curr_BW(ihots)/sum(curr_BW);
    end
    
    [Ie_zE] = Ie_M_stream_4_aurora(h_atm,mag_ze,E(1:iEmax),mu_lims,mu_scatterings,Ie_oftGB,p_e_q,ne,Te,...
                                   nO,O_levels,XsO(:,1:iEmax),@O_e_2nd_dist,@phase_fcn_O,...
                                   nN2,N2_levels,XsN2(:,1:iEmax),@N2_e_2nd_dist,@phase_fcn_N2,...
                                   nO2,O2_levels,XsO2(:,1:iEmax),@O2_e_2nd_dist,@phase_fcn_O2);
    
    savefile = fullfile(savedir,sprintf('IeSteady_4-%02d.mat',i1));
    save(savefile,'Ie_zE','E','mu_lims','h_atm','Ie_oftGB','mu_scatterings','-v7.3')
  end
  
  
end % Run_section(1)

