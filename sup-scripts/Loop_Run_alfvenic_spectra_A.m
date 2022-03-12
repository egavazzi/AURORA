
Science_grade_not_demo = 1;

if Science_grade_not_demo
  n_dirs = 1 + 180*4
  stepE = 1
else % Then demo with coarser energy-grid
  n_dirs = 1 + 180*4;
  stepE = 30;
end
print_figs = 0;
plot_figs = 1;
%% Time intervall
t = linspace(0,1.5,151);

Isteady = 1e12;
I_flick = 1e12;

if ~exist('setup_completed','var') ||  setup_completed ~= 1
  setup4etrptd10streamsB2 % setup4etrptdms
end
disp('push any key')
pause
% pars for AWacc:   z0,  dtheta@s, n_t   BeamWeights    scE  Dir-name
%                  (m)         (m)  #                    C   char
par_list4_AWacc = {2400e3,     10,  8, 2*pi*B_W(1:5),    1, 'AWA_2400-10-1';
                   2400e3,     10,  8, 2*pi*B_W(1:5),    2, 'AWA_2400-10-2';
                   2400e3,     10,  8, 2*pi*B_W(1:5),    4, 'AWA_2400-10-4';
                   3000e3,     10,  8, 2*pi*B_W(1:5),    1, 'AWA_3000-10-1';
                   3000e3,     10,  8, 2*pi*B_W(1:5),    2, 'AWA_3000-10-2';
                   3000e3,     10,  8, 2*pi*B_W(1:5),    4, 'AWA_3000-10-4';
                   2400e3,     20,  8, 2*pi*B_W(1:5),    1, 'AWA_2400-20-1';
                   2400e3,     20,  8, 2*pi*B_W(1:5),    2, 'AWA_2400-20-2';
                   2400e3,     20,  8, 2*pi*B_W(1:5),    4, 'AWA_2400-20-4';
                   3000e3,     20,  8, 2*pi*B_W(1:5),    1, 'AWA_3000-20-1';
                   3000e3,     20,  8, 2*pi*B_W(1:5),    2, 'AWA_3000-20-2';
                   3000e3,     20,  8, 2*pi*B_W(1:5),    4, 'AWA_3000-20-4'};
                   t_off = 1.65; % was: 0.65
I00 = zeros(numel(h_atm)*(numel(mu_lims)-1),numel(E));

load('Ie_Lailas_snug.mat','Ec','tC','Ie_Bpar')
tC0 = tC;
Ec0 = Ec;
method = 'linear'; % Or 'nearest'
%% Found a suitable Gaussian flux at higher energy
%E0dE = 8*[1,0.1]*1e3;
%EfdEf = [E0dE(1)+E0dE(2)*1,E0dE(2)*4.374];

%% High-Energy steady-state precipitation 3000 km, 0-15 deg, 8.8 keV
%
%  This section calculates the rise-characteristics of
%  electron-fluxes in the ionosphere for 8.8 keV mono-energetic
%  field-aligned precipitation (0-15 degrees pitch-angle).
% %% Etrp_Flickering_HET4_9stream.m:savedir = 'Onset_H8p8keV

dE = 11.65; % TODO: Fix this properly, whatever that means
theta0 = (0:89)*pi/180;
mu0 = cos(theta0);
for i_pars = 2:size(par_list4_AWacc,1) %Run_section(1)
  disp('===========================================================================')
  disp(par_list4_AWacc(i_pars,:))
  disp('---------------------------------------------------------------------------')
  dtheta_at_source = par_list4_AWacc{i_pars,2}*pi/180;  % 2390e3 3390e3; % Source altitude at 3000 km
  z_0      = par_list4_AWacc{i_pars,1}; % 0, 5, 10
  n_loop   = par_list4_AWacc{i_pars,3}; % 5, 7;
  curr_BW  = par_list4_AWacc{i_pars,4}; % B_W;
  scE      = par_list4_AWacc{i_pars,5}  % ;
  savedir  = par_list4_AWacc{i_pars,6}  % 'Onset_H8keV';
  curr_par = par_list4_AWacc(i_pars,:)
  
  Ec = Ec0*scE;
  Eref = max(Ec);
  %curr_BW = curr_BW.*(1/3).^(0:(numel(curr_BW)-1))
  [dEmax,iEmax] = min(abs(E-Eref))
  %% 0 Hertz flickering 8 keV
  I0 = 0*squeeze(I00(:,1:iEmax));
  [SUCCESS,MESSAGE,MESSAGEID] = mkdir(savedir);
  save(fullfile(savedir,'neutral_atm.mat'),'Te','h_atm','nN2','nO','nO2','ne')
  save(fullfile(savedir,'curr_par.mat'),'curr_par')
  
  % E0 = maxE;% TODO - fix this thing
  E0 = max(Ec);% TODO - fix this thing
  if 1
    [B,lB,lon,lat,alt,x_B,y_B,z_B] = make_BofL([],[],[],h_atm(end),z_0);
    for i1 = 1:90,
      [t_delay,mu_source] = dt_of_z_mirror(v_of_E(Eref),lB(1),lB(end),cos((i1-1)*pi/180),B,lB);
      dt0(i1) = t_delay;                                                                     
      MU_at_source(i1) = mu_source;                                                          
    end
    dt00 = min(dt0);
  end
  
  
  fprintf('%d: %s\n',i1,datestr(now,'HH:MM:SS'))
  
  p_e_q = zeros(length(h_atm),length(E(1:iEmax)));
  % disp([E0dE(1,1),E0dE(1,2)])
  
  for iBeam = numel(c_o_mu):-1:1
    Ie_oftE{iBeam} = @(t,E) zeros(numel(E),numel(t));
  end
  i_Hots = 1:numel(curr_BW);
  for ihots = i_Hots,
    % Ie_oftGB{ihots}=@(t,Ei) Ie_ParMod(t,z_source,z_0,Ei,abs(real(c_o_mu(ihots))),t_off,11.65,scE);
    % Ie_oftGB{ihots} = @(t,Ei) Ie_tmuE_par(t,abs(real(c_o_mu(ihots))),Ei,...
    %                                       dt0,MU_at_source,E0,mu0,t_off,dE,scE);
    Ie_oftE{ihots} = @(ti,Ei) Ie_tmuE_par(tC,Ec,cos((0:89)*pi/180),Ie_Bpar,...
                                          ti,Ei,...
                                          -c_o_mu(ihots),...
                                          Eref,dt0,MU_at_source,dt00,cos(dtheta_at_source),median(dE),method)*curr_BW(ihots);    
  end
  for i1 = 1:n_loop,
    
    t_run = t(1:21)+t(21)*(i1-1);
    [Ie_ztE] = Ie_Mstream_tz_2_aurora(h_atm,mag_ze,E(1:iEmax),mu_lims,mu_scatterings,...
                                      t_run,...
                                      I0,Ie_oftE,p_e_q,ne,Te,...
                                      [],...
                                      nO,O_levels,XsO(:,1:iEmax),@O_e_2nd_dist,@phase_fcn_O,...
                                      nN2,N2_levels,XsN2(:,1:iEmax),@N2_e_2nd_dist,@phase_fcn_N2,...
                                      nO2,O2_levels,XsO2(:,1:iEmax),@O2_e_2nd_dist,@phase_fcn_O2);
    
    I0 = squeeze(Ie_ztE(:,end,:));
    savefile = fullfile(savedir,sprintf('IeFlickering-%02d.mat',i1));
    save(savefile,'Ie_ztE','E','t_run','mu_lims','h_atm','Ie_oftE','I0','mu_scatterings','-v7.3')
  end
  
end % Run_section(1)


for i_pars = 1:size(par_list4_AWacc,1) %Run_section(1)
  disp('===========================================================================')
  disp(par_list4_AWacc(i_pars,:))
  disp('---------------------------------------------------------------------------')
  dtheta_at_source = par_list4_AWacc{i_pars,2}*pi/180;  % 2390e3 3390e3; % Source altitude at 3000 km
  z_0      = par_list4_AWacc{i_pars,1}; % 0, 5, 10
  n_loop   = par_list4_AWacc{i_pars,3}; % 5, 7;
  curr_BW  = par_list4_AWacc{i_pars,4}; % B_W;
  scE      = par_list4_AWacc{i_pars,5}  % ;
  savedir  = par_list4_AWacc{i_pars,6}  % 'Onset_H8keV';
  curr_par = par_list4_AWacc(i_pars,:)
  
  Ec = Ec0*scE;
  Eref = max(Ec);
  %curr_BW = curr_BW.*(1/3).^(0:(numel(curr_BW)-1))
  [dEmax,iEmax] = min(abs(E-Eref))
  %% 0 Hertz flickering 8 keV
  I0 = 0*squeeze(I00(:,1:iEmax));
  [SUCCESS,MESSAGE,MESSAGEID] = mkdir(savedir);
  save(fullfile(savedir,'neutral_atm.mat'),'Te','h_atm','nN2','nO','nO2','ne')
  save(fullfile(savedir,'curr_par.mat'),'curr_par')
  
  % E0 = maxE;% TODO - fix this thing
  E0 = max(Ec);% TODO - fix this thing
  if 1
    [B,lB,lon,lat,alt,x_B,y_B,z_B] = make_BofL([],[],[],h_atm(end),z_0);
    for i1 = 1:90,
      [t_delay,mu_source] = dt_of_z_mirror(v_of_E(Eref),lB(1),lB(end),cos((i1-1)*pi/180),B,lB);
      dt0(i1) = t_delay;                                                                     
      MU_at_source(i1) = mu_source;                                                          
    end
    dt00 = min(dt0);
  end
  
  
    
  p_e_q = zeros(length(h_atm),length(E(1:iEmax)));
  % disp([E0dE(1,1),E0dE(1,2)])
  
  for iBeam = numel(c_o_mu):-1:1
    Ie_oftE{iBeam} = @(t,E) zeros(numel(E),numel(t));
  end
  i_Hots = 1:numel(curr_BW);
  for ihots = i_Hots,
    % Ie_oftGB{ihots}=@(t,Ei) Ie_ParMod(t,z_source,z_0,Ei,abs(real(c_o_mu(ihots))),t_off,11.65,scE);
    % Ie_oftGB{ihots} = @(t,Ei) Ie_tmuE_par(t,abs(real(c_o_mu(ihots))),Ei,...
    %                                       dt0,MU_at_source,E0,mu0,t_off,dE,scE);
    Ie_oftE{ihots} = @(ti,Ei) Ie_tmuE_par(tC,Ec,cos((0:89)*pi/180),Ie_Bpar,...
                                          ti,Ei,...
                                          -c_o_mu(ihots),...
                                          Eref,dt0,MU_at_source,dt00,cos(dtheta_at_source),median(dE),method)*curr_BW(ihots);    
  end
  
  for ihots = 1:10,
    for iE = 1:iEmax,   
      Ie_test(:,iE,ihots) = Ie_oftE{ihots}(t,E(iE));  
    end
  end
  figure
  if 1                                                                                                                                        
    subplot(4,1,1)
    pcolor(t,E(1:iEmax),log10(squeeze(Ie_test(:,1:iEmax,1)))'),shading flat,caxis([-4 0]+max(caxis)),set(gca,'yscale','log'),colorbar
    axis([0 1.5 30 E(iEmax)])
    subplot(4,1,2)
    pcolor(t,180-acos(c_o_mu)*180/pi,log10(squeeze(Ie_test(:,iEmax-20,:)))'),shading flat,caxis([-4 0]+max(caxis)),colorbar
    title(E(iEmax-20))
    subplot(4,1,3)
    pcolor(t,180-acos(c_o_mu)*180/pi,log10(squeeze(Ie_test(:,iEmax-55,:)))'),shading flat,caxis([-4 0]+max(caxis)),colorbar
    title(E(iEmax-55))
    subplot(4,1,4)
    pcolor(t,180-acos(c_o_mu)*180/pi,log10(squeeze(Ie_test(:,iEmax-70,:)))'),shading flat,caxis([-4 0]+max(caxis)),colorbar
    title(E(iEmax-70))
  end
  
end
