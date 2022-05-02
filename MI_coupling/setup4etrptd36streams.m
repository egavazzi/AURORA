%% Set-up script for Examples in time-dependent electron transport

%% Figure-generation
plot_figs = 0;
print_figs = 0;
if print_figs
  mkdir('Figures')
end

%% Loading atmosphere model
%
% The electron transport needs the density profiles of the major
% thermospheric constituents, to model the collisions (ellastic and in-
% ellastc) excitations and ionizations. Here we use the MSIS-00 atmosphere
% for the date: 
load msis20051008_3.dat

h_atm = (74:1:610)'*1e3;h_atm = (100:1:610)'*1e3;
dz = @(n) 150 + 150/200*(0:(n-1))' +1.2*exp(((0:(n-1))-150)/17)';
z_atm = 100e3 + cumsum(dz(308)) - dz(1);
numel(z_atm)/numel(h_atm);
h_atm = z_atm(:); % make sure it's a column vector

OPS.atmosphere = interp1(msis20051008_3(:,1),msis20051008_3,h_atm/1e3,'pchip');
OPS.atmosphere(:,2:4) = OPS.atmosphere(:,2:4)*1e6; % into m^-3
OPS.atmosphere(:,5) = OPS.atmosphere(:,5)*1e3;     % into kg/m^-3
% Below is an attempt to make a gradual transition from the dilute
% uppermost thermosphere to mathematical vacuum, this might be useful to
% ensure that the gradient of the upward fluxes automatically becomes
% mathematically zero at the top to be in agreement with the mathematical
% boundary condition for the upward fluxes. This is one possible cause of
% peculiar numerical behaviour of the solution
nO = OPS.atmosphere(:,2);
nN2 = OPS.atmosphere(:,3);
nO2 = OPS.atmosphere(:,4);
nO(end-2:end)  = 0;
nO2(end-2:end) = 0;
nN2(end-2:end) = 0;
%nO(end-7:end-3)  = (erf((2:-1:-2)'/2)+1)/2.*nO(end-7:end-3);
%nO2(end-7:end-3) = (erf((2:-1:-2)'/2)+1)/2.*nO2(end-7:end-3);
%nN2(end-7:end-3) = (erf((2:-1:-2)'/2)+1)/2.*nN2(end-7:end-3);
nO(end-5:end-3)  = (erf((1:-1:-1)'/2)+1)/2.*nO(end-5:end-3);
nO2(end-5:end-3) = (erf((1:-1:-1)'/2)+1)/2.*nO2(end-5:end-3);
nN2(end-5:end-3) = (erf((1:-1:-1)'/2)+1)/2.*nN2(end-5:end-3);

%%
% Plot the atmosphere
if plot_figs
  figure_atm = figure;
  subplot(1,2,1)
  phMSIS = semilogx([nN2,nO2,nO],h_atm/1e3,'-','linewidth',2);
  axis([1e12 1e20 85 620])
  ylabel('height (km)')
  xlabel('[m^{-3}]')
  title('MSIS neutral density')
  legend(phMSIS,'N_2','O_2','O')
  set(phMSIS(2),'linestyle','--')
  set(phMSIS(3),'linestyle','-.')
end
%% 1.5 Load n_e(z) as estimated from EISCAT observations
%  In addition to N_2, O_2 and O we need an electron density profile, to
%  account for the energy degradation of primary electrons caused by the
%  accumulated effect of millions of very weak Coulomb collisions, ideally
%  ne(z) should be gotten from ISR-data, but for testing this is OK
load iri20051008.dat
Iri20051008 = iri20051008(1:end,:);
ne = interp1(Iri20051008(:,1),Iri20051008(:,2),h_atm/1e3,'pchip')*10;
ne(h_atm<82e3) = 1;
%% Plot electron density profile
if plot_figs
  figure(figure_atm)
  subplot(1,2,2)
  semilogx(ne,h_atm/1e3,'linewidth',2)
  title('electron density')
  xlabel('[m^{-3}]')
  axis([1e8 1e13 85 520])
  set(gca,'xtick',10.^(6:15))
  if print_figs
    filename = fullfile('Figures','MSIS-IRI-atm-01.eps');
    [status, message] = backup1000filesversions(filename);
    print('-depsc2','-painters',filename)
  end
end
% Te = IRI!
Te = interp1(Iri20051008(:,1),Iri20051008(:,6),h_atm/1e3,'pchip');
Te(Te==-1) = 350;


%% Magnetic inclination
mag_ze = 13*pi/180; % This is the approximate magnetic inclination at EISCAT

%% Excitation thresholds
% And finally we load the energy levels of the excited states of
% the major neutral species, this is obviously necessary to model the
% energy-losses in inellastic collisions leading to excitation of discrete
% levels.
load N2_levels.dat
load O2_levels.dat
load O_levels.dat
% Plot energy levels
if plot_figs
  figure
  subplot(1,3,1)
  plot([0 1],O_levels(:,1)*[1 1],'-','linewidth',2)
  axis([0 1 0 25])
  set(gca,'xtick',[])
  ylabel('Excitation energy (eV)')
  title('O')
  subplot(1,3,2)
  plot([0 1],O2_levels(:,1)*[1 1],'-','linewidth',2)
  axis([0 1 0 25])
  set(gca,'xtick',[])
  title('O_2')
  subplot(1,3,3)
  plot([0 1],N2_levels(:,1)*[1 1],'-','linewidth',2)
  axis([0 1 0 25])
  set(gca,'xtick',[])
  title('N_2')
  if print_figs
    filename = fullfile('Figures','ON2O2-Energy-levels-01.eps');
    [status, message] = backup1000filesversions(filename);
    print('-depsc2','-painters',filename)
    % print('-depsc2','-painters','ON2O2-Energy-levels-01.eps')
  end
end
%% Energy grids
% For the transport-functions t work properly it seems as if the
% energy-grid should not have bins wider than the lowest ionization
% threshold. 
%
% This energy-grid below work well.
dEfcn = @(X,DEinitial,DEfinal,C,X0) DEinitial+(1+tanh(C*(X-X0)))/2*DEfinal;
E = cumsum(dEfcn(0:1600,0.15,11.5,0.05,80))+1.9;
E = E(1:stepE:end);
dE = diff(E);dE = [dE,dE(end)];
if plot_figs
  figure
  plot(E,dE,'r.-')
  xlabel('energy [eV]')
  ylabel('energy [eV]')
  title('Energy variation of energy bin size')
  if print_figs
    filename = fullfile('Figures','Energy-grid-01.eps');
    [status, message] = backup1000filesversions(filename);
    print('-depsc2','-painters',filename)
    % print('-depsc2','-painters','Energy-grid-01.eps')
  end
end
%% Collision cross sections
%  These collision cross-sections are taken from a series of
%  compilation-papers by Itikawa and colleagues, from there the
%  cross-sections have been digitized. (Itikawa et al. 1986 (N2), Itikawa
%  and Ichimura 1990 (O), and Itikawa et al 1989 (O2))
[XsO,xs_fcnO] =   get_all_xs('O',E+dE/2);
[XsO2,xs_fcnO2] = get_all_xs('O2',E+dE/2);
[XsN2,xs_fcnN2] = get_all_xs('N2',E+dE/2);
% %% Back-scattering-ratios
% %  For two-stream calculations the angular part of the scattering is simply
% %  an energy-dependent fraction of electron that escapes a collision in the
% %  opposite hemispherical direction after a collision. These ratios are
% %  taken from Solomon et al. 1988
% 
% bscrO = back_scattering_42stream('O',E+dE/2);
% bscrO2 = back_scattering_42stream('O2',E+dE/2);
% bscrN2 = back_scattering_42stream('N2',E+dE/2);
%% Name/-/label of excited levels.
fid = fopen('N2_levels.name','r');
qwe = textscan(fid,'%s');
fclose(fid);
iN = 1;
for i1 = 1:numel(qwe{1})
  if qwe{1}{i1}(1) ~= '%'
    N2names{iN} = qwe{1}{i1};
    iN = iN+1;
  end
end
% N2names = qwe{1};
fid = fopen('O2_levels.name','r');
qwe = textscan(fid,'%s');
fclose(fid);
iN = 1;
for i1 = 1:numel(qwe{1})
  if qwe{1}{i1}(1) ~= '%'
    O2names{iN} = qwe{1}{i1};
    iN = iN+1;
  end
end
%O2names = qwe{1};
fid = fopen('O_levels.name','r');
qwe = textscan(fid,'%s');
fclose(fid);
iN = 1;
for i1 = 1:numel(qwe{1})
  if qwe{1}{i1}(1) ~= '%'
    Onames{iN} = qwe{1}{i1};
    iN = iN+1;
  end
end
% Onames = qwe{1};

%% Plot cross-sections and back-scattering ratios.
%  For illustration and illumination, note that at energies larger than a
%  few 10s of eV the ionization-cross-section dominates over the other
%  in-ellastic cross-sections
if plot_figs
  figure('position',[306, 415, 1746, 559])
  subplot(1,3,1)
  phO = loglog(E,XsO,'linewidth',2);
  axis([E(1) E(end) 1e-25 1e-18])
  axis([2.0539 8000 1e-24 3e-19])
  set(phO(1),'color','k')
  set(phO(2:7),'linestyle','-')
  colormap(jet(64))
  cmlines(phO(2:7))
  set(phO(8:end),'linestyle','--')
  cmpC = inferno(64);
  colormap(cmpC(1:60,:))
  cmlines(phO(8:end))
  % legend(phO,Onames) % TODO: strrep fine_%d1_%d2 with fine(%d1-%d2)
  for i1 = numel(phO):-1:1
    Onames4l{i1} = sprintf('%d',i1);
  end
  legend(phO,Onames4l) % TODO: strrep fine_%d1_%d2 with fine(%d1-%d2)
  set(gca,'xtick',[1 10 100 1000 1e4])
  title('\sigma O','fontsize',16)
  ylabel('[m^{-2}]','fontsize',16)
  xlabel('[eV]','fontsize',16)
  set(gca,'fontsize',14)
  
  subplot(1,3,2)
  phN2 = loglog(E,XsN2,'linewidth',2);
  axis([E(1) E(end) 1e-25 1e-18])
  axis([2.0539 8000 1e-24 3e-19])
  cmlines(phN2)
  set(phN2(2:3:end),'linestyle','--')
  set(phN2(3:3:end),'linestyle','-.')
  set(phN2(end),'color','k')
  for i1 = numel(phN2):-1:1
    N2names4l{i1} = sprintf('%d',i1);
  end
  %legend(phN2,N2names4l)
  % legend(phN2,N2names)
  set(phN2(1),'color','k')
  set(phN2(2:12),'linestyle','-','linewidth',1) 
  % colormap(summer)
  colormap(py_D_cmap)
  cmlines(phN2(2:12))     
  set(phN2(end-5:end),'linestyle','--')         
  set(phN2(13:end-6),'linestyle','-')  
  cmpC = inferno(64);
  colormap(cmpC(1:60,:))
  cmlines(phN2(end-5:end))
  colormap(jet(64))
  cmlines(phN2(13:end-6)) 
  title('\sigma N_2','fontsize',16)
  set(gca,'xtick',[1 10 100 1000 1e4])
  xlabel('[eV]','fontsize',16)
  %legHdl = gridLegend(phN2,2,N2names4l);
  set(gca,'fontsize',14)
  [legend_h,object_h,plot_h,text_str] = legendflex(phN2,N2names4l,'ncol',2);  
  set(legend_h,'position',get(legend_h,'position')+[40 20 0 0])
  subplot(1,3,3)
  phO2 = loglog(E,XsO2([1,9,8,7,6,5,4,3,2,10:15],:),'linewidth',2);
  axis([E(1) E(end) 1e-25 1e-18])
  axis([2.0539 8000 1e-24 3e-19])
  set(phO2(2),'linewidth',1)
  colormap(py_D_cmap)
  set(phO2(2),'color',[0 0.5 0])  
  colormap(jet(64))
  cmlines(phO2(3:9))
  set(phO2(10:end),'linestyle','--')
  cmpC = inferno(64);
  colormap(cmpC(1:60,:))
  cmlines(phO2(10:end))
  set(phO2(1),'color','k')
  for i1 = numel(phO2):-1:1
    O2names4l{i1} = sprintf('%d',i1);
  end
  legend(phO2,O2names4l)
  title('\sigma O_2','fontsize',16)
  set(gca,'xtick',[1 10 100 1000 1e4])
  xlabel('[eV]','fontsize',15)
  set(gca,'fontsize',14)
  if print_figs
    filename = fullfile('Figures','ON2O2-cross-sections-01.eps');
    [status, message] = backup1000filesversions(filename);
    print('-depsc2','-painters',filename)
    % print('-depsc2','-painters','ON2O2-cross-sections-01.eps')
  end
  
  figure % N2-only
  set(gcf,'position',[31, 348, 1209, 746])
  phN2 = loglog(E,XsN2,'linewidth',2);
  axis([E(1) E(end) 1e-25 1e-18])
  axis([2.0539 8000 1e-24 3e-19])
  cmlines(phN2)
  set(phN2(2:3:end),'linestyle','--')
  set(phN2(3:3:end),'linestyle','-.')
  set(phN2(end),'color','k')
  for i1 = numel(phN2):-1:1
    N2names4l{i1} = sprintf('%d',i1);
  end
  %legend(phN2,N2names4l)
  % legend(phN2,N2names)
  set(phN2(1),'color','k')
  set(phN2(2:12),'linestyle','-','linewidth',1) 
  % colormap(summer)
  colormap(py_D_cmap)
  cmlines(phN2(2:12))     
  set(phN2(end-5:end),'linestyle','--')         
  set(phN2(13:end-6),'linestyle','-')  
  cmpC = inferno(64);
  colormap(cmpC(1:60,:))
  cmlines(phN2(end-5:end))
  colormap(jet(64))
  cmlines(phN2(13:end-6)) 
  title('\sigma N_2','fontsize',16)
  set(gca,'xtick',[1 10 100 1000 1e4])
  xlabel('[eV]','fontsize',16)
  %legHdl = gridLegend(phN2,2,N2names4l);
  set(gca,'fontsize',14)
  
  legHdl = gridLegend(phN2,2,N2names);
  if print_figs
    filename = fullfile('Figures','N2-cross-sections-01.eps');
    [status, message] = backup1000filesversions(filename);
    print('-depsc2','-painters',filename)
    % print('-depsc2','-painters','N2-cross-sections-01.eps')
  end
  figure % O2-only
  set(gcf,'position',[31, 348, 1209, 746])
  phO2 = loglog(E,XsO2([1,9,8,7,6,5,4,3,2,10:15],:),'linewidth',2);
  axis([E(1) E(end) 1e-25 1e-18])
  axis([2.0539 8000 1e-24 3e-19])
  set(phO2(2),'linewidth',1)
  colormap(py_D_cmap)
  set(phO2(2),'color',[0 0.5 0])  
  colormap(jet(64))
  cmlines(phO2(3:9))
  set(phO2(10:end),'linestyle','--')
  cmpC = inferno(64);
  colormap(cmpC(1:60,:))
  cmlines(phO2(10:end))
  set(phO2(1),'color','k')
  for i1 = numel(phO2):-1:1
    O2names4l{i1} = sprintf('%d',i1);
  end
  legend(phO2,O2names([1,9,8,7,6,5,4,3,2,10:15]))
  title('\sigma O_2','fontsize',16)
  set(gca,'xtick',[1 10 100 1000 1e4])
  xlabel('[eV]','fontsize',15)
  set(gca,'fontsize',14)

  if print_figs
    filename = fullfile('Figures','O2-cross-sections-01.eps');
    [status, message] = backup1000filesversions(filename);
    print('-depsc2','-painters',filename)
    % print('-depsc2','-painters','N2-cross-sections-01.eps')
  end
  figure % O - only
%  [~,idxOlevels] = sort(O_levels(:,1));
idxOlevels = 1:size(O_levels,1);
%   O_levels = O_levels(idxOlevels,:);
%   Onames = Onames(idxOlevels);
  set(gcf,'position',[31, 348, 1209, 746])
  phO = loglog(E,XsO(idxOlevels,:),'linewidth',2);
  set(gca,'fontsize',12)
  cmlines(phO)
  title('\sigma O')
  xlabel('[eV]','fontsize',16)
  ylabel('[m^{-2}]','fontsize',16)
  axis([E(1) E(end) 1e-25 1e-18])
  set(gca,'xtick',[1 10 100 1000 1e4])
  set(phO(end-3:end),'linestyle','--')         
  set(phO(1:end-4),'linestyle','-','linewidth',2)  
  cmlines(phO(1:end-4)) 
  cmlines(phO(end-3:end)) 
  legHdl = gridLegend(phO,1,Onames(idxOlevels));
  if print_figs
    filename = fullfile('Figures','O-cross-sections-01.eps');
    [status, message] = backup1000filesversions(filename);
    print('-depsc2','-painters',filename)
    % print('-depsc2','-painters','N2-cross-sections-01.eps')
  end
  
end
%% Pre-calculations of cascadning electron-spectra for ionizations
S2ndO = O_e_2nd_dist(E,E(end),O_levels(end,1),'c',AURORA_root_directory);
S2ndO2 = O2_e_2nd_dist(E,E(end),O2_levels(end,1),'c',AURORA_root_directory);
S2ndN2 = N2_e_2nd_dist(E,E(end),N2_levels(end,1),'c',AURORA_root_directory);


%% 10-stream mu-limits:
theta_lims2do = 180:-5:0;
%% Calculate the beam-to-beam scattering and weigthings
%  This can be done in the function Ie_Mstream_tz_2_aurora, as was done
%  above where the mu_scatterings argument was set to a place-holder 1, or
%  before as done here. This is necessary to get the proper weights for
%  calculating the fluxes in different beams to more for example an
%  isotropic flux
[Pmu2mup,theta2beamW,BeamW,mu_lims] = e_scattering_result_finder(theta_lims2do,AURORA_root_directory);
mu_scatterings = {Pmu2mup,theta2beamW,BeamW};
c_o_mu = mu_avg(mu_lims);
% B_W is created for easify generation of functions for the
% downward fluxes and upward flux gradients at the upper boundary,
% i.e. the downward fluxes are the primary electron spectra, and
% since the boundary condition for the upward fluxes are that their
% gradient is zero it is nice to make a variable that guarantees that
B_W = BeamW;
B_W(c_o_mu>0) = 0;
if plot_figs
  figure('position',[20, 495, 1220, 458])
  for i1 = numel(BeamW):-1:1,
    theta_str{i1} = sprintf('%3.1f^\\circ - %3.1f^\\circ',...
                            180-180/pi*acos(mu_lims(i1)),...
                            180-180/pi*acos(mu_lims(i1+1)));
  end
  spp = [2*ones(10,1),5*ones(10,1),[1:5,(10:-1:6)]'];
  for i1 = 1:10,
    subplot(spp(i1,1),spp(i1,2),spp(i1,3)) 
    imagesc(0:180,180:-1:0,Pmu2mup(:,:,i1)),caxis([0 1])          
    title(theta_str{i1},'fontsize',14)
    if i1 == 8
      xlabel('scattering-angle','fontsize',14)
    elseif i1 == 5 || i1 == 6
      colorbar_labeled('')
    elseif i1 == 1 || i1 == 10
      ylabel('from-angle','fontsize',14)
    end
    set(gca,'xtick',0:30:180,'ytick',0:30:180,'tickdir','out','box','off')
  end
  colormap(inferno)
  
  if print_figs
    filename = fullfile('Figures','Scattering-Matrices-02.eps');
    [status, message] = backup1000filesversions(filename);
    print('-depsc2','-painters',filename)
    % print('-depsc2','-painters','Scattering-Matrices-02.eps')
  end
  figure
  imagesc(180:-1:0,1:10,theta2beamW)
  title('Beam-weigths','fontsize',15)
  xlabel('pitch-angle (degrees)','fontsize',15)
  ylabel('pitch-angle-stream','fontsize',15)
  colorbar
  if print_figs
    filename = fullfile('Figures','Beam-Weights-02.eps');
    [status, message] = backup1000filesversions(filename);
    print('-depsc2','-painters',filename)
    % print('-depsc2','-painters','Beam-Weights-02.eps')
  end
end

%%

I_stable = 1e12;

setup_completed = 1;
