%% 12-stream

stepE = 1;
setup4etrptd10streams

t = 0:0.005:0.05;

Isteady = 1e12;
I_flick = 1e12;


pause(0.3)
theta_lims2do = fliplr([0 10 25 45 65 80 90 100 115 135 155 170 180])

[Pmu2mup,theta2beamW,BeamW,mu_lims] = e_scattering_result_finder(theta_lims2do,...
                                                  AURORA_root_directory);
c_o_mu = mu_avg(mu_lims);

% load e_s_b_18streams.mat theta_lims mu_lims Pmu2mup theta2beamW BeamW c_o_mu
mu_scatterings = {Pmu2mup,theta2beamW,BeamW};
B_W = BeamW.*(c_o_mu<0)

I00 = zeros(numel(h_atm)*(numel(mu_lims)-1),numel(E));

% Energy grid for the low-energy tail
E_LET = E;
% Initially we use many runs without the LET, so simply set a
% zero-flux in all pitch-angle-streams:
for i1 = (numel(mu_lims)-1):-1:1,
  Ie_LET{i1} = 0*E_LET;
end


par_list4G_fa = {30e3, 0.0001,   6,        B_W(1),  [1e3 100], 'G_fa_st12-3000-1keV-0-10'
                 30e3, 0.0001,   6,      B_W(1:6),  [1e3 100], 'G_fa_st12-3000-1keV-0-90'};


%%
%  This section calculates the rise-characteristics of
%  electron-fluxes in the ionosphere for mono-energetic
%  Precipitation, with energy and pitch-angle distribution as
%  parameterized in PAR_LIST4G_FA


for i_pars = 1:size(par_list4G_fa,1) %Run_section(1)
  disp('===========================================================================')
  disp(par_list4G_fa(i_pars,:))
  disp('---------------------------------------------------------------------------')
  z_source = par_list4G_fa{i_pars,1};  % 2390e3 3390e3; % Source altitude at 3000 km
  f_f      = par_list4G_fa{i_pars,2}; % 0, 5, 10
  n_loop   = par_list4G_fa{i_pars,3}; % 5, 7;
  curr_BW  = par_list4G_fa{i_pars,4} % B_W;
  E0dE     = par_list4G_fa{i_pars,5}  % ;
  savedir  = par_list4G_fa{i_pars,6}  % 'Onset_H8keV';
  curr_par = par_list4G_fa(i_pars,:)
  
  h_over_zmax = z_source - h_atm(end);
  %curr_BW = curr_BW.*(1/3).^(0:(numel(curr_BW)-1))
  iEmax = find(    cumsum(exp(-(E-E0dE(1)).^2/E0dE(2)^2)) < ...
               0.9999*sum(exp(-(E-E0dE(1)).^2/E0dE(2)^2)),1,'last')
  %% 0 Hertz flickering 3 keV
  I0 = 0*squeeze(I00(:,1:iEmax));
  [SUCCESS,MESSAGE,MESSAGEID] = mkdir(savedir);
  save(fullfile(savedir,'neutral_atm.mat'),'Te','h_atm','nN2','nO','nO2','ne')
  save(fullfile(savedir,'curr_par.mat'),'curr_par')
  t_run = t;
  
  for i1 = 1:n_loop,
    
    fprintf('%d: %s\n',i1,datestr(now,'HH:MM:SS'))
    
    p_e_q = zeros(length(h_atm),length(E(1:iEmax)));
    % disp([E0dE(1,1),E0dE(1,2)])
    
    for iBeam = numel(c_o_mu):-1:1
      Ie_oftGB{iBeam} = @(t,E) Ie_flickeringGauss_fluxS(t,E,...
                                                       0*Isteady,E0dE(1,1),E0dE(1,2),...
                                                       Ie_LET{iBeam},E_LET,...
                                                       f_f*pi*2,...
                                                       500e3,E(end),0)*B_W(iBeam).*ones(size(t));
    end
    i_Hots = 1:numel(curr_BW);
    for ihots = i_Hots,
      Ie_oftGB{ihots} = @(t,Ei) Ie_flickeringGauss_fluxS(t,Ei,...
                                                    Isteady,E0dE(1,1),E0dE(1,2),...
                                                    Ie_LET{ihots},E_LET,...
                                                    f_f*pi*2,...
                                                    h_over_zmax,E(end),...
                                                    0*I_flick,E0dE(1),E0dE(2),...
                                                    abs(real(c_o_mu(1))),1)*curr_BW(ihots)/sum(curr_BW).*ones(size(t));
    end
    [Ie_ztE] = Ie_Mstream_tz_2_aurora(h_atm,mag_ze,E(1:iEmax),mu_lims,mu_scatterings,...
                                      t_run,...
                                      I0,Ie_oftGB,p_e_q,ne,Te,...
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



%% 18-stream

theta_lims2do = [180, 175, 167, 156, 143, 127, 114, 103,  95,  90, ...
                 85,  77,  66,  53,  37,  24,  13,   5,   0]; 
[Pmu2mup,theta2beamW,BeamW,mu_lims] = e_scattering_result_finder(theta_lims2do,...
                                                  AURORA_root_directory);
c_o_mu = mu_avg(mu_lims);

% load e_s_b_18streams.mat theta_lims mu_lims Pmu2mup theta2beamW BeamW c_o_mu
mu_scatterings = {Pmu2mup,theta2beamW,BeamW};
B_W = BeamW.*(c_o_mu<0)

par_list4G_fa = {30e3, 0.0001,   6,      B_W(1:9),  [1e3 100], 'G_fa_st18-3000-1keV-0-90'};

I00 = zeros(numel(h_atm)*(numel(mu_lims)-1),numel(E));

% Energy grid for the low-energy tail
E_LET = E;
% Initially we use many runs without the LET, so simply set a
% zero-flux in all pitch-angle-streams:
for i1 = (numel(mu_lims)-1):-1:1,
  Ie_LET{i1} = 0*E_LET;
end


%%
%  This section calculates the rise-characteristics of
%  electron-fluxes in the ionosphere for mono-energetic
%  Precipitation, with energy and pitch-angle distribution as
%  parameterized in PAR_LIST4G_FA

for i_pars = 1:size(par_list4G_fa,1) %Run_section(1)
  disp('===========================================================================')
  disp(par_list4G_fa(i_pars,:))
  disp('---------------------------------------------------------------------------')
  z_source = par_list4G_fa{i_pars,1};  % 2390e3 3390e3; % Source altitude at 3000 km
  f_f      = par_list4G_fa{i_pars,2}; % 0, 5, 10
  n_loop   = par_list4G_fa{i_pars,3}; % 5, 7;
  curr_BW  = par_list4G_fa{i_pars,4} % B_W;
  E0dE     = par_list4G_fa{i_pars,5}  % ;
  savedir  = par_list4G_fa{i_pars,6}  % 'Onset_H8keV';
  curr_par = par_list4G_fa(i_pars,:)
  
  h_over_zmax = z_source - h_atm(end);
  %curr_BW = curr_BW.*(1/3).^(0:(numel(curr_BW)-1))
  iEmax = find(    cumsum(exp(-(E-E0dE(1)).^2/E0dE(2)^2)) < ...
               0.9999*sum(exp(-(E-E0dE(1)).^2/E0dE(2)^2)),1,'last')
  %% 0 Hertz flickering 3 keV
  I0 = 0*squeeze(I00(:,1:iEmax));
  [SUCCESS,MESSAGE,MESSAGEID] = mkdir(savedir);
  save(fullfile(savedir,'neutral_atm.mat'),'Te','h_atm','nN2','nO','nO2','ne')
  save(fullfile(savedir,'curr_par.mat'),'curr_par')
  t_run = t;
  
  for i1 = 1:n_loop,
    
    fprintf('%d: %s\n',i1,datestr(now,'HH:MM:SS'))
    
    p_e_q = zeros(length(h_atm),length(E(1:iEmax)));
    % disp([E0dE(1,1),E0dE(1,2)])
    
    for iBeam = numel(c_o_mu):-1:1
      Ie_oftGB{iBeam} = @(t,E) Ie_flickeringGauss_fluxS(t,E,...
                                                       0*Isteady,E0dE(1,1),E0dE(1,2),...
                                                       Ie_LET{iBeam},E_LET,...
                                                       f_f*pi*2,...
                                                       500e3,E(end),0)*B_W(iBeam).*ones(size(t));
    end
    i_Hots = 1:numel(curr_BW);
    for ihots = i_Hots,
      Ie_oftGB{ihots} = @(t,Ei) Ie_flickeringGauss_fluxS(t,Ei,...
                                                    Isteady,E0dE(1,1),E0dE(1,2),...
                                                    Ie_LET{ihots},E_LET,...
                                                    f_f*pi*2,...
                                                    h_over_zmax,E(end),...
                                                    0*I_flick,E0dE(1),E0dE(2),...
                                                    abs(real(c_o_mu(1))),1)*curr_BW(ihots)/sum(curr_BW).*ones(size(t));
    end
    [Ie_ztE] = Ie_Mstream_tz_2_aurora(h_atm,mag_ze,E(1:iEmax),mu_lims,mu_scatterings,...
                                      t_run,...
                                      I0,Ie_oftGB,p_e_q,ne,Te,...
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


%% 18-stream-2

theta_lims2do = [180, 170, 160, 150, 140, 130, 120, 110, 100,  90, ...
                 80,  70,  60,  50,  40,  30,  20,  10,   0];
[Pmu2mup,theta2beamW,BeamW,mu_lims] = e_scattering_result_finder(theta_lims2do,...
                                                  AURORA_root_directory);
c_o_mu = mu_avg(mu_lims);

% load e_s_b_19_streams_20191122-004141.mat theta_lims mu_lims Pmu2mup theta2beamW BeamW c_o_mu
mu_scatterings = {Pmu2mup,theta2beamW,BeamW};
B_W = BeamW.*(c_o_mu<0)

par_list4G_fa = {30e3, 0.0001,   6,        B_W(1),  [1e3 100], 'G_fa_st19-3000-1keV-0-10'
                 30e3, 0.0001,   6,      B_W(1:9),  [1e3 100], 'G_fa_st19-3000-1keV-0-90'};


I00 = zeros(numel(h_atm)*(numel(mu_lims)-1),numel(E));


% Energy grid for the low-energy tail
E_LET = E;
% Initially we use many runs without the LET, so simply set a
% zero-flux in all pitch-angle-streams:
for i1 = (numel(mu_lims)-1):-1:1,
  Ie_LET{i1} = 0*E_LET;
end

%%
%  This section calculates the rise-characteristics of
%  electron-fluxes in the ionosphere for mono-energetic
%  Precipitation, with energy and pitch-angle distribution as
%  parameterized in PAR_LIST4G_FA

for i_pars = 1:size(par_list4G_fa,1) %Run_section(1)
  disp('===========================================================================')
  disp(par_list4G_fa(i_pars,:))
  disp('---------------------------------------------------------------------------')
  z_source = par_list4G_fa{i_pars,1};  % 2390e3 3390e3; % Source altitude at 3000 km
  f_f      = par_list4G_fa{i_pars,2}; % 0, 5, 10
  n_loop   = par_list4G_fa{i_pars,3}; % 5, 7;
  curr_BW  = par_list4G_fa{i_pars,4} % B_W;
  E0dE     = par_list4G_fa{i_pars,5}  % ;
  savedir  = par_list4G_fa{i_pars,6}  % 'Onset_H8keV';
  curr_par = par_list4G_fa(i_pars,:)
  
  h_over_zmax = z_source - h_atm(end);
  %curr_BW = curr_BW.*(1/3).^(0:(numel(curr_BW)-1))
  iEmax = find(    cumsum(exp(-(E-E0dE(1)).^2/E0dE(2)^2)) < ...
               0.9999*sum(exp(-(E-E0dE(1)).^2/E0dE(2)^2)),1,'last')
  %% 0 Hertz flickering 3 keV
  I0 = 0*squeeze(I00(:,1:iEmax));
  [SUCCESS,MESSAGE,MESSAGEID] = mkdir(savedir);
  save(fullfile(savedir,'neutral_atm.mat'),'Te','h_atm','nN2','nO','nO2','ne')
  save(fullfile(savedir,'curr_par.mat'),'curr_par')
  t_run = t;
  
  for i1 = 1:n_loop,
    
    fprintf('%d: %s\n',i1,datestr(now,'HH:MM:SS'))
    
    p_e_q = zeros(length(h_atm),length(E(1:iEmax)));
    % disp([E0dE(1,1),E0dE(1,2)])
    
    for iBeam = numel(c_o_mu):-1:1
      Ie_oftGB{iBeam} = @(t,E) Ie_flickeringGauss_fluxS(t,E,...
                                                       0*Isteady,E0dE(1,1),E0dE(1,2),...
                                                       Ie_LET{iBeam},E_LET,...
                                                       f_f*pi*2,...
                                                       500e3,E(end),0)*B_W(iBeam).*ones(size(t));
    end
    i_Hots = 1:numel(curr_BW);
    for ihots = i_Hots,
      Ie_oftGB{ihots} = @(t,Ei) Ie_flickeringGauss_fluxS(t,Ei,...
                                                    Isteady,E0dE(1,1),E0dE(1,2),...
                                                    Ie_LET{ihots},E_LET,...
                                                    f_f*pi*2,...
                                                    h_over_zmax,E(end),...
                                                    0*I_flick,E0dE(1),E0dE(2),...
                                                    abs(real(c_o_mu(1))),1)*curr_BW(ihots)/sum(curr_BW).*ones(size(t));
    end
    [Ie_ztE] = Ie_Mstream_tz_2_aurora(h_atm,mag_ze,E(1:iEmax),mu_lims,mu_scatterings,...
                                      t_run,...
                                      I0,Ie_oftGB,p_e_q,ne,Te,...
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


%% 20-stream


theta_lims2do = [180, 176, 169, 159, 148, 135, 122, 111, 101,  94, ...
                 90,  86,  79,  69,  58,  45,  32,  21,  11,   4,   0];
[Pmu2mup,theta2beamW,BeamW,mu_lims] = e_scattering_result_finder(theta_lims2do,...
                                                  AURORA_root_directory);
c_o_mu = mu_avg(mu_lims);
% $$$ load e_s_b_20streams.mat theta_lims mu_lims Pmu2mup theta2beamW BeamW c_o_mu
mu_scatterings = {Pmu2mup,theta2beamW,BeamW};
B_W = BeamW.*(c_o_mu<0)
I00 = zeros(numel(h_atm)*(numel(mu_lims)-1),numel(E));


par_list4G_fa = {30e3, 0.0001,   6,      B_W(1:10),  [1e3 100], 'G_fa_st20-3000-1keV-0-90'};

% Energy grid for the low-energy tail
E_LET = E;
% Initially we use many runs without the LET, so simply set a
% zero-flux in all pitch-angle-streams:
for i1 = (numel(mu_lims)-1):-1:1,
  Ie_LET{i1} = 0*E_LET;
end

%%
%  This section calculates the rise-characteristics of
%  electron-fluxes in the ionosphere for mono-energetic
%  Precipitation, with energy and pitch-angle distribution as
%  parameterized in PAR_LIST4G_FA

for i_pars = 1:size(par_list4G_fa,1) %Run_section(1)
  disp('===========================================================================')
  disp(par_list4G_fa(i_pars,:))
  disp('---------------------------------------------------------------------------')
  z_source = par_list4G_fa{i_pars,1};  % 2390e3 3390e3; % Source altitude at 3000 km
  f_f      = par_list4G_fa{i_pars,2}; % 0, 5, 10
  n_loop   = par_list4G_fa{i_pars,3}; % 5, 7;
  curr_BW  = par_list4G_fa{i_pars,4} % B_W;
  E0dE     = par_list4G_fa{i_pars,5}  % ;
  savedir  = par_list4G_fa{i_pars,6}  % 'Onset_H8keV';
  curr_par = par_list4G_fa(i_pars,:)
  
  h_over_zmax = z_source - h_atm(end);
  %curr_BW = curr_BW.*(1/3).^(0:(numel(curr_BW)-1))
  iEmax = find(    cumsum(exp(-(E-E0dE(1)).^2/E0dE(2)^2)) < ...
               0.9999*sum(exp(-(E-E0dE(1)).^2/E0dE(2)^2)),1,'last')
  %% 0 Hertz flickering 3 keV
  I0 = 0*squeeze(I00(:,1:iEmax));
  [SUCCESS,MESSAGE,MESSAGEID] = mkdir(savedir);
  save(fullfile(savedir,'neutral_atm.mat'),'Te','h_atm','nN2','nO','nO2','ne')
  save(fullfile(savedir,'curr_par.mat'),'curr_par')
  t_run = t;
  
  for i1 = 1:n_loop,
    
    fprintf('%d: %s\n',i1,datestr(now,'HH:MM:SS'))
    
    p_e_q = zeros(length(h_atm),length(E(1:iEmax)));
    % disp([E0dE(1,1),E0dE(1,2)])
    
    for iBeam = numel(c_o_mu):-1:1
      Ie_oftGB{iBeam} = @(t,E) Ie_flickeringGauss_fluxS(t,E,...
                                                       0*Isteady,E0dE(1,1),E0dE(1,2),...
                                                       Ie_LET{iBeam},E_LET,...
                                                       f_f*pi*2,...
                                                       500e3,E(end),0)*B_W(iBeam).*ones(size(t));
    end
    i_Hots = 1:numel(curr_BW);
    for ihots = i_Hots,
      Ie_oftGB{ihots} = @(t,Ei) Ie_flickeringGauss_fluxS(t,Ei,...
                                                    Isteady,E0dE(1,1),E0dE(1,2),...
                                                    Ie_LET{ihots},E_LET,...
                                                    f_f*pi*2,...
                                                    h_over_zmax,E(end),...
                                                    0*I_flick,E0dE(1),E0dE(2),...
                                                    abs(real(c_o_mu(1))),1)*curr_BW(ihots)/sum(curr_BW).*ones(size(t));
    end
    [Ie_ztE] = Ie_Mstream_tz_2_aurora(h_atm,mag_ze,E(1:iEmax),mu_lims,mu_scatterings,...
                                      t_run,...
                                      I0,Ie_oftGB,p_e_q,ne,Te,...
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
