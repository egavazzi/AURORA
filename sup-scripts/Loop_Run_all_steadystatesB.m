% $$$ function Ie_zE = Loop_Run_all_steadystatesB()
%% Set-up stuff
% Here we set variables for the running of this script
Science_grade_not_demo = 1; % To run with an OK energy-grid
plot_figs = 1;              % To plot intermediate figures during
                            % set-up, to disable set to zero
print_figs = 0;             % To print the intermediate figures
                            % during set-up, change to one

if Science_grade_not_demo
  n_dirs = 1 + 180*4;
  stepE = 1;
else % Then demo with coarser energy-grid
  n_dirs = 1 + 180*4;
  stepE = 30;
end

%% Time intervall


%Isteady = 1e12;

if ~exist('setup_completed','var') ||  setup_completed ~= 1
  setup4etrptd10streamsB
end
par_list4Ie = {3000e3, 0.0,   1,         1,  [5e2 100], 'G_Beam1-0.5keV';
               3000e3, 0.0,   1,         1,  [1e3 100], 'G_Beam1-1.0keV';
               3000e3, 0.0,   1,         1,  [2e3 100], 'G_Beam1-2.0keV';
               3000e3, 0.0,   1,         1,  [4e3 100], 'G_Beam1-4.0keV';
               3000e3, 0.0,   1,         1,  [6e3 100], 'G_Beam1-6.0keV';
               3000e3, 0.0,   1,       1:5,  [5e2 100], 'G_IS015-0.5keV';
               3000e3, 0.0,   1,       1:5,  [1e3 100], 'G_ISO15-1.0keV';
               3000e3, 0.0,   1,       1:5,  [2e3 100], 'G_ISO15-2.0keV';
               3000e3, 0.0,   1,       1:5,  [4e3 100], 'G_ISO15-4.0keV';
               3000e3, 0.0,   1,       1:5,  [6e3 100], 'G_ISO15-6.0keV'};

par_list4Ie = { B,               1,  [1e3 100], 'Bm_Beam1-1keV-1';
                B,               2,  [1e3 100], 'Bm_Beam2-1keV-1';
                B,               3,  [1e3 100], 'Bm_Beam3-1keV-1';
                B,               4,  [1e3 100], 'Bm_Beam4-1keV-1';
                B,               5,  [1e3 100], 'Bm_Beam5-1keV-1';
                B,             1:5,  [1e3 100], 'Bm_ISO15-1keV-1'
                ones(size(B)),   1,  [1e3 100], 'nB_Beam1-1keV-1';
                ones(size(B)),   2,  [1e3 100], 'nB_Beam2-1keV-1';
                ones(size(B)),   3,  [1e3 100], 'nB_Beam3-1keV-1';
                ones(size(B)),   4,  [1e3 100], 'nB_Beam4-1keV-1';
                ones(size(B)),   5,  [1e3 100], 'nB_Beam5-1keV-1';
                ones(size(B)), 1:5,  [1e3 100], 'nB_ISO15-1keV-1'
  
              };

               
for i_pars = 1:size(par_list4Ie,1) %Run_section(1)
  disp('===========================================================================')
  disp(par_list4Ie(i_pars,:))
  disp('---------------------------------------------------------------------------')
  B_curr   = par_list4Ie{i_pars,1}; % Magnetic field;
  idx_BW   = par_list4Ie{i_pars,2}; % pitch-angle-stream-index(es)
  E0dE     = par_list4Ie{i_pars,3}  % Peak energy and width of Gaussian
  savedir  = par_list4Ie{i_pars,4}  % Directory for saving results to
  
  curr_par = par_list4Ie(i_pars,:);

  
  iEmax = find( cumsum(exp(-(E-E0dE(1)).^2/E0dE(2)^2)) < ...
                0.9999*sum(exp(-(E-E0dE(1)).^2/E0dE(2)^2)),1,'last');
  
  [SUCCESS,MESSAGE,MESSAGEID] = mkdir(savedir);
  save(fullfile(savedir,'neutral_atm.mat'),'Te','h_atm','nN2','nO','nO2','ne')
  save(fullfile(savedir,'curr_par.mat'),'curr_par')
  
  fprintf('%d: %s\n',i_pars,datestr(now,'HH:MM:SS'))
  
  p_e_q = zeros(length(h_atm),length(E(1:iEmax)));
  
  for iBeam = numel(c_o_mu):-1:1
    Ie_oftGB{iBeam} = @(E) zeros(size(E));
  end
  i_Hots = idx_BW;
  for i_curr_Beam = i_Hots,
    disp(i_curr_Beam)
    Ie_oftGB{i_curr_Beam} = @(Ei) 1e12*exp(-(Ei-E0dE(1,1)).^2/E0dE(1,2))*B_W(i_curr_Beam)/sum(B_W(i_Hots));
  end
  
  [Ie_zE] = Ie_M_stream_4_auroraB(h_atm,mag_ze,B_curr,E(1:iEmax),mu_lims,mu_scatterings,Ie_oftGB,p_e_q,ne,Te,...
                                  nO,O_levels,XsO(:,1:iEmax),@O_e_2nd_dist,@phase_fcn_O,...
                                  nN2,N2_levels,XsN2(:,1:iEmax),@N2_e_2nd_dist,@phase_fcn_N2,...
                                  nO2,O2_levels,XsO2(:,1:iEmax),@O2_e_2nd_dist,@phase_fcn_O2);
  
  savefile = fullfile(savedir,sprintf('IeSteady_1.mat'));
  save(savefile,'Ie_zE','E','mu_lims','h_atm','Ie_oftGB','mu_scatterings','curr_par','-v7.3')
  Ie_zE_all{i_pars} = Ie_zE;
end



%% Results and analysis 1 electron-fluxes
%% 1, electron fluxes as functions of altitude and energy
% pitch-angle-stream by pitch-angle-stream

fignames = {'stream #1 B-field-mirror-force',
            'stream #2 B-field-mirror-force',
            'stream #3 B-field-mirror-force',
            'stream #4 B-field-mirror-force',
            'stream #5 B-field-mirror-force',
            'streams #1-5 B-field-mirror-force',
            'stream #1 no mirror-force',
            'stream #2 no mirror-force',
            'stream #3 no mirror-force',
            'stream #4 no mirror-force',
            'stream #5 no mirror-force',
            'streams #1-5 no mirror-force'};
for i1 = 1:12
  figure('name',fignames{i1},'position',[104 534 1186 435])
  plot_IezE_3DEz(h_atm,E(1:size(Ie_zE_all{i1},2)),Ie_zE_all{i1},dE(1:size(Ie_zE_all{i1},2)),BeamW,[6 12],spp, theta_str);   
end
%% 2 electron fluxes at 1 keV
% as function of pitch-angle-stream-number and altitude
[dE1000,iE1000] = min(abs(E-1e3));
figure('name','Ie(E=1keV,z,\mu),B-mirror')
for i1 = 1:6,
  Ie1000 = reshape(Ie_zE_all{i1}(:,iE1000),numel(h_atm),10)./repmat(BeamW*2*pi*dE(iE1000),size(h_atm));
  subplot(2,3,i1)
  pcolor(0.5:10.5,h_atm/1e3,log10(Ie1000(:,[1:end,end]))),shading flat,caxis([-6 0]+max(caxis))
  cx(i1+6,:) = caxis;
  set(gca,'xtick',1:10,'tickdir','out','box','off')
  if i1 == 3 || i1 == 6
    colorbar_labeled('e^-/m^2/s/ster/eV','log')
  else
    colorbar_labeled('','log')
  end
  if i1 == 1 || i1 == 4
    ylabel('height (km)')
  else
    set(gca,'yticklabel','')
  end
  if i1 > 3
    xlabel('pitch-angle-stream #')
  end
end

figure('name','Ie(E=1keV,z,\mu),no-mirror')
for i1 = 1:6,
  Ie1000 = reshape(Ie_zE_all{i1+6}(:,iE1000),numel(h_atm),10)./repmat(BeamW*2*pi*dE(iE1000),size(h_atm));
  subplot(2,3,i1)
  pcolor(0.5:10.5,h_atm/1e3,log10(Ie1000(:,[1:end,end]))),shading flat,caxis([-6 0]+max(caxis))
  cx(i1+6,:) = caxis;
  set(gca,'xtick',1:10,'tickdir','out','box','off')
  if i1 == 3 || i1 == 6
    colorbar_labeled('e^-/m^2/s/ster/eV','log')
  else
    colorbar_labeled('','log')
  end
  if i1 == 1 || i1 == 4
    ylabel('height (km)')
  else
    set(gca,'yticklabel','')
  end
  if i1 > 3
    xlabel('pitch-angle-stream #')
  end
end
%% 3 Pitch-angle spectra at highest altitude
figure
subplot(2,3,1)
plot_IezE_2DEtheta(h_atm/1e3,E,dE,Ie_zE_all{1},BeamW,mu_lims,numel(h_atm)-1,[]);
subplot(2,3,2)
plot_IezE_2DEtheta(h_atm/1e3,E,dE,Ie_zE_all{4},BeamW,mu_lims,numel(h_atm)-1,[]); 
subplot(2,3,3)
plot_IezE_2DEtheta(h_atm/1e3,E,dE,Ie_zE_all{6},BeamW,mu_lims,numel(h_atm)-1,[]);
subplot(2,3,4)
plot_IezE_2DEtheta(h_atm/1e3,E,dE,Ie_zE_all{1+6},BeamW,mu_lims,numel(h_atm)-1,[]);
subplot(2,3,5)
plot_IezE_2DEtheta(h_atm/1e3,E,dE,Ie_zE_all{4+6},BeamW,mu_lims,numel(h_atm)-1,[]);
subplot(2,3,6)
plot_IezE_2DEtheta(h_atm/1e3,E,dE,Ie_zE_all{6+6},BeamW,mu_lims,numel(h_atm)-1,[]);

figure('name','Ie(z=600 km,E,\mu),B-mirror')
E_flux_up = zeros(12,1);
E_flux_down = zeros(12,1);

for i1 = 1:6
  Ie_top = Ie_zE_all{i1}(numel(h_atm):numel(h_atm):end,:)./(BeamW'* ...
                                                  dE(1:size(Ie_zE_all{1},2))*2*pi);
  Ie_topraw = Ie_zE_all{i1}(numel(h_atm):numel(h_atm):end,:);
  Ie_up(i1,:) = sum(Ie_topraw(6:10,:));
  Ie_down(i1,:) = sum(Ie_topraw(1:5,:));
  subplot(2,3,i1)
  pcolor(E(1:size(Ie_zE_all{1},2)),0.5:10.5,log10(Ie_top([1:end,end],:)))
  shading flat
  caxis([-6 0]+max(caxis))
  set(gca,'ytick',1:10,'yticklabel',theta_str)
  if i1 == 3 || i1 == 6
    colorbar_labeled('e^-/m^2/s/ster/eV','log')
  else
    colorbar_labeled('','log')
  end
  if i1 == 1 || i1 == 4
    set(gca,'ytick',1:10,'yticklabel',theta_str)
  else
    set(gca,'ytick',1:10,'yticklabel','')
  end
  if i1 > 3
    xlabel('Energy')
  end
  
  for i2 = 6:10
    E_flux_up(i1) = E_flux_up(i1) + trapz(E(1:size(Ie_zE_all{1},2)),E(1:size(Ie_zE_all{1},2)).*Ie_topraw(i2,:));
  end
  for i2 = 1:5
    E_flux_down(i1) = E_flux_down(i1) + trapz(E(1:size(Ie_zE_all{1},2)),E(1:size(Ie_zE_all{1},2)).*Ie_topraw(i2,:));
  end
end
figure('name','Ie(z=600 km,E,\mu),no-B-mirror')
for i1 = 1:6
  Ie_top = Ie_zE_all{i1+6}(numel(h_atm):numel(h_atm):end,:)./(BeamW'* ...
                                                  dE(1:size(Ie_zE_all{1},2))*2*pi);
  Ie_topraw = Ie_zE_all{i1+6}(numel(h_atm):numel(h_atm):end,:);
  Ie_up(i1+6,:) = sum(Ie_topraw(6:10,:));
  Ie_down(i1+6,:) = sum(Ie_topraw(1:5,:));
  subplot(2,3,i1)
  pcolor(E(1:size(Ie_zE_all{1},2)),0.5:10.5,log10(Ie_top([1:end,end],:)))
  shading flat
  caxis([-6 0]+max(caxis))
  set(gca,'ytick',1:10,'yticklabel',theta_str)
  if i1 == 3 || i1 == 6
    colorbar_labeled('e^-/m^2/s/ster/eV','log')
  else
    colorbar_labeled('','log')
  end
  if i1 == 1 || i1 == 4
    set(gca,'ytick',1:10,'yticklabel',theta_str)
  else
    set(gca,'ytick',1:10,'yticklabel','')
  end
  if i1 > 3
    xlabel('Energy')
  end
  
  for i2 = 6:10
    E_flux_up(i1+6) = E_flux_up(i1+6) + trapz(E(1:size(Ie_zE_all{1},2)),E(1:size(Ie_zE_all{1},2)).*Ie_topraw(i2,:));
  end
  for i2 = 1:5
    E_flux_down(i1+6) = E_flux_down(i1+6) + trapz(E(1:size(Ie_zE_all{1},2)),E(1:size(Ie_zE_all{1},2)).*Ie_topraw(i2,:));
  end
end



%% Excitation and ionization-rates

emXS4278 = exc_4278(E);
emXS6730 = exc_6730_N2(E);
emXS8446_O = exc_8446_O(E);
emXS8446_O2 = exc_8446_O2(E);
emX7774_O = exc_7774_O(E);
emX7774_O2 = exc_7774_O2(E);
excXS_O1D = exc_O1D(E);
excXS_O1S = exc_O1S(E);

dE = diff(E);
dE = dE([1:end,end]);

XsO  = get_all_xs('O',E+dE/2);
XsO2 = get_all_xs('O2',E+dE/2);
XsN2 = get_all_xs('N2',E+dE/2);

load N2_levels.dat
load O2_levels.dat
load O_levels.dat

XsOi = O_levels(:,2)'*XsO;
XsO2i = O2_levels(:,2)'*XsO2;
XsN2i = N2_levels(:,2)'*XsN2;
XsN2A3 = XsN2(13,:);


for i1 = 1:numel(Ie_zE_all)
  
  szIzE = size(Ie_zE_all{i1});
  Q4278(:,i1) = exc_z_of_Ie_zE(h_atm,E(1:szIzE(2)),Ie_zE_all{i1},nN2,emXS4278(1:szIzE(2)));
  
  Q6730(:,i1) = exc_z_of_Ie_zE(h_atm,E(1:szIzE(2)),Ie_zE_all{i1},nN2,emXS6730(1:szIzE(2)));
  
  Q8446(:,i1) = exc_z_of_Ie_zE(h_atm,E(1:szIzE(2)),Ie_zE_all{i1},nO,emXS8446_O(1:szIzE(2))) + ...
          exc_z_of_Ie_zE(h_atm,E(1:szIzE(2)),Ie_zE_all{i1},nO2,emXS8446_O2(1:szIzE(2)));
  Q7774(:,i1) = exc_z_of_Ie_zE(h_atm,E(1:szIzE(2)),Ie_zE_all{i1},nO,emX7774_O(1:szIzE(2))) + ...
          exc_z_of_Ie_zE(h_atm,E(1:szIzE(2)),Ie_zE_all{i1},nO2,emX7774_O2(1:szIzE(2)));
  
  Q8446_O(:,i1) = exc_z_of_Ie_zE(h_atm,E(1:szIzE(2)),Ie_zE_all{i1},nO,emXS8446_O(1:szIzE(2)));
  Q7774_O(:,i1) = exc_z_of_Ie_zE(h_atm,E(1:szIzE(2)),Ie_zE_all{i1},nO,emX7774_O(1:szIzE(2)));
  
  Q8446_O2(:,i1) = exc_z_of_Ie_zE(h_atm,E(1:szIzE(2)),Ie_zE_all{i1},nO2,emXS8446_O2(1:szIzE(2)));
  Q7774_O2(:,i1) = exc_z_of_Ie_zE(h_atm,E(1:szIzE(2)),Ie_zE_all{i1},nO2,emX7774_O2(1:szIzE(2)));
  
  QO1D(:,i1) = exc_z_of_Ie_zE(h_atm,E(1:szIzE(2)),Ie_zE_all{i1},nO,excXS_O1D(1:szIzE(2)));
  QO1S(:,i1) = exc_z_of_Ie_zE(h_atm,E(1:szIzE(2)),Ie_zE_all{i1},nO,excXS_O1S(1:szIzE(2)));
  
  QN2i(:,i1) = exc_z_of_Ie_zE(h_atm,E(1:szIzE(2)),Ie_zE_all{i1},nN2,XsN2i(1:szIzE(2)));
  QO2i(:,i1) = exc_z_of_Ie_zE(h_atm,E(1:szIzE(2)),Ie_zE_all{i1},nO2,XsO2i(1:szIzE(2)));
  QOi(:,i1) = exc_z_of_Ie_zE(h_atm,E(1:szIzE(2)),Ie_zE_all{i1},nO,XsOi(1:szIzE(2)));
  
  QN2A3(:,i1) = exc_z_of_Ie_zE(h_atm,E(1:szIzE(2)),Ie_zE_all{i1},nN2,XsN2A3(1:szIzE(2)));
  
end

