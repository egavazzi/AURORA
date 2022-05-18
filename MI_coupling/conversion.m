% This script converts the distribution function extracted from the
% MI-model and expressed as a function of vz and mu_mag into a distribution
% function expressed as a function of E and mu_pitch that can be used as
% input in a ionosphere response model.

%% Import data
inputb6;
%load fzvzmu3950000.mat
load Bfield.mat
%% Calculate flux ([#e/m2/s]) for specie 1 (magnetospheric e-)
zz = 600;     % index for the altitude at top of ionosphere

% Extracting fzvzmu_in %still importing data
cd ./f_vlasov_raw/
dDir = dir;
f_in = zeros(36,100,100);
load /mnt/data/etienne/ketchup/MI_coupling/l7/1.27e7_21.5to120s_results/outp/fzvzmu3950000.mat
f_in(1,:,:) = fzvzmustruct(1).f(particle(1).vz >= 0,:,zz);
for i1 = 3:numel(dDir)
  load(dDir(i1).name);
  f_in(i1-1,:,:) = fzvzmustruct(1).f(particle(1).vz >= 0,:,zz);    % taking only downward going e-
  disp(['file ',num2str(i1-1),'/',num2str(size(f_in,1)),' done!'])
end
cd ../
%%
tic
index_specie = 1;
m = particle(index_specie).mass;   
B = B(end);   % magnetic field at top of ionosphere

% Extraction of the (vz,mu_mag)-grid
vz_middle_bin = particle(index_specie).vz(particle(index_specie).vz >= 0);      % taking only downward going e-
vz_grid = particle(index_specie).vzcorn(particle(index_specie).vzcorn >= 0);    % taking only downward going e-
dvz = particle(index_specie).dvz;
mu_mag_middle_bin = particle(index_specie).mu;
mu_mag_grid = particle(index_specie).mucorn;
dmu_mag = particle(index_specie).dmu;

% Extrapolate f over a finer (vz,mu_mag) grid;
HMR_VZ = 20;   % HOW MUCH DO YOU WANT TO REFINE VZ
HMR_MU = 20;   % HOW MUCH DO YOU WANT TO REFINE MU
% Refine vz-grid
F = griddedInterpolant(1:101,vz_grid);
vz_grid_finer = F(1:(1/HMR_VZ):101);
dvz_finer = (1/HMR_VZ) * dvz ;
vz_middle_bin_finer = vz_grid_finer(2:end) - 0.5 * dvz_finer;
% Refine mu_mag-grid
G = griddedInterpolant(1:101,mu_mag_grid);
mu_mag_grid_finer = G(1:(1/HMR_MU):101);
dmu_mag_finer = (1/HMR_MU) .* dmu_mag;
dmu_mag_finer = repelem(dmu_mag_finer,HMR_MU);
mu_mag_middle_bin_finer = mu_mag_grid_finer(2:end) - 0.5 .* dmu_mag_finer;


% Make E-grid
dEfcn = @(X,DEinitial,DEfinal,C,X0) DEinitial+(1+tanh(C*(X-X0)))/2*DEfinal;     % from setup4etrpd10streams.m
E_grid = cumsum(dEfcn(0:700,0.15,11.5,0.05,80))+1.9;
E_middle_bin = E_grid;
dE = diff(E_grid);dE = [dE,dE(end)];
% Make mu_pitch-grid
theta_lims2do = [180:-10:90];
mu_pitch_grid = cosd(theta_lims2do);
BeamW = [];
for iMu = (numel(mu_pitch_grid)-1):-1:1
  BeamW(iMu) = 2*pi * abs(integral(@(pa) sin(pa),acos(mu_pitch_grid(iMu)),acos(mu_pitch_grid(iMu+1))));
end
mu_pitch_middle_bin = mu_avg(mu_pitch_grid);


Ie_1_out = zeros(size(f_in,1),length(E_middle_bin),length(mu_pitch_middle_bin));
for i_t = 1:size(f_in,1)
  % Extraction of the distr funct f(v,mu_mag)
  fzvzmu_in = f_in(i_t,:,:);
  fzvzmu_in = squeeze(fzvzmu_in);
  % Refine f
  fzvzmu_in_finer = repelem(fzvzmu_in,HMR_VZ,HMR_MU);
  
  Ie_1 = zeros(numel(E_middle_bin),numel(mu_pitch_middle_bin));
  v = zeros(length(vz_middle_bin),length(mu_mag_middle_bin));
  for ii = 1:length(vz_middle_bin)
    for jj = 1:length(mu_mag_middle_bin)
      v(ii,jj) = sqrt(vz_middle_bin(ii).^2 + 2*B/m .* mu_mag_middle_bin(jj));
    end
  end
  v = repelem(v,HMR_VZ,HMR_MU);

  for ii = 1:length(vz_middle_bin_finer)
    for jj = 1:length(mu_mag_middle_bin_finer)
      % get coordinate point in (E,mu_pitch)-grid
      E = 0.5 * m * (vz_middle_bin_finer(ii).^2 + 2*B/m .* mu_mag_middle_bin_finer(jj)) ./ 1.6e-19;
      mu_pitch = - cos(atan(sqrt(2*B/m * mu_mag_middle_bin_finer(jj)) / vz_middle_bin_finer (ii)));
      [~,index_pitch] = min(abs(mu_pitch_middle_bin - mu_pitch));
      [~,index_energy] = min(abs(E_middle_bin - E));
      % convert distribution function and compute flux
      Ie_1(index_energy,index_pitch) = Ie_1(index_energy,index_pitch) + ...
                                  v(ii,jj) * fzvzmu_in_finer(ii,jj) * ...
                                  dvz_finer * dmu_mag_finer(jj);
    end
  end
%   Ie_total_1 = sum(((v .* fzvzmu_in_finer) * dmu_mag_finer.') * dvz_finer)    % [#e/m2/s]
%   E_total_1 = sum(E_middle_bin * (Ie_1 ))                                     % [eV/m2/s]
  Ie_1_out(i_t,:,:) = Ie_1;
  disp(['Conversion ',num2str(i_t),'/',num2str(size(f_in,1)),' done!'])
end
toc
%% Calculate flux ([#e/m2/s]) for specie 3 (ionospheric e-)
tic
index_specie = 3;
m = particle(index_specie).mass;
zz = 300;      % index for the altitude at top of ionosphere
B = B(end);   % magnetic field at top of ionosphere

% Extraction of the (vz,mu_mag)-grid
vz_middle_bin = particle(index_specie).vz(particle(index_specie).vz <= 0);      % taking only downward going e-
vz_grid = particle(index_specie).vzcorn(particle(index_specie).vzcorn <= 0);    % taking only downward going e-
dvz = particle(index_specie).dvz;
mu_mag_middle_bin = particle(index_specie).mu;
mu_mag_grid = particle(index_specie).mucorn;
dmu_mag = particle(index_specie).dmu;
% Extraction of the distr funct f(v,mu_mag)
fzvzmu_in = fzvzmustruct(index_specie).f(particle(index_specie).vz <= 0,:,zz);  % taking only downward going e-

% Extrapolate f over a finer (vz,mu_mag) grid;
HMR_VZ = 20;   % HOW MUCH DO YOU WANT TO REFINE VZ
HMR_MU = 20;   % HOW MUCH DO YOU WANT TO REFINE MU
% Refine vz-grid
F = griddedInterpolant(1:101,vz_grid);
vz_grid_finer = F(1:(1/HMR_VZ):101);
dvz_finer = (1/HMR_VZ) * dvz ;
vz_middle_bin_finer = vz_grid_finer(2:end) - 0.5 * dvz_finer;
% Refine mu_mag-grid
G = griddedInterpolant(1:101,mu_mag_grid);
mu_mag_grid_finer = G(1:(1/HMR_MU):101);
dmu_mag_finer = (1/HMR_MU) .* dmu_mag;
dmu_mag_finer = repelem(dmu_mag_finer,HMR_MU);
% mu_mag_middle_bin_finer = mu_mag_middle_bin;
mu_mag_middle_bin_finer = mu_mag_grid_finer(2:end) - 0.5 .* dmu_mag_finer;
% Refine f
fzvzmu_in_finer = repelem(fzvzmu_in,HMR_VZ,HMR_MU);

% Make E-grid
dEfcn = @(X,DEinitial,DEfinal,C,X0) DEinitial+(1+tanh(C*(X-X0)))/2*DEfinal;     % from setup4etrpd10streams.m
E_grid = cumsum(dEfcn(0:1000,0.15,11.5,0.05,80))+1.9;
E_middle_bin = E_grid;
dE = diff(E_grid);dE = [dE,dE(end)];
% Make mu_pitch-grid
theta_lims2do = [180:-5:90];
mu_pitch_grid = cosd(theta_lims2do);
BeamW = [];
for iMu = (numel(mu_pitch_grid)-1):-1:1
  BeamW(iMu) = 2*pi * abs(integral(@(pa) sin(pa),acos(mu_pitch_grid(iMu)),acos(mu_pitch_grid(iMu+1))));
end
mu_pitch_middle_bin = mu_avg(mu_pitch_grid);

Ie_3 = zeros(numel(E_middle_bin),numel(mu_pitch_middle_bin));

v = zeros(length(vz_middle_bin),length(mu_mag_middle_bin));
for ii = 1:length(vz_middle_bin)
  for jj = 1:length(mu_mag_middle_bin)
    v(ii,jj) = sqrt(vz_middle_bin(ii).^2 + 2*B/m .* mu_mag_middle_bin(jj));
  end
end
v = repelem(v,HMR_VZ,HMR_MU);

for ii = 1:length(vz_middle_bin_finer)
  for jj = 1:length(mu_mag_middle_bin_finer)
    % get coordinate point in (E,mu_pitch)-grid
    E = 0.5 * m * (vz_middle_bin_finer(ii).^2 + 2*B/m .* mu_mag_middle_bin_finer(jj)) ./ 1.6e-19;
    mu_pitch = - cos(atan(sqrt(2*B/m * mu_mag_middle_bin_finer(jj)) / vz_middle_bin_finer (ii)));
    [~,index_pitch] = min(abs(mu_pitch_middle_bin - mu_pitch));
    [~,index_energy] = min(abs(E_middle_bin - E));
    % convert distribution function and compute flux (/eV/ster)
    Ie_3(index_energy,index_pitch) = Ie_3(index_energy,index_pitch) + ...
                                v(ii,jj) * fzvzmu_in_finer(ii,jj) * ...
                                dvz_finer * dmu_mag_finer(jj);
  end
end
Ie_total_3 = sum(((v .* fzvzmu_in_finer) * dmu_mag_finer.') * dvz_finer)    % [#e/m2/s]
E_total_3 = sum(E_middle_bin * (Ie_3))                                      % [eV/m2/s]
toc
%% In which bins is the bulk of the energy flux
E1 = E_middle_bin.' .* Ie_1;
R = sum(sum(E1(log10(E1)<2)))./sum(sum(E1)) % log10(E)>1 gives around 100% of the energy. Means we can reduce size of E-grid.
%% plot Ie
Ie_plot = E1;%Ie_1;
% Ie_plot = Ie_3;
% Ie_plot = Ie_1 + Ie_3; 

figure()
c1 = jet(64); c2 = jet(256); c3 = jet(1024);
c = [c1(1:8,:);c2(32:64,:);c3(257:1024,:)];
colormap(c)
pp = log10(Ie_plot);
ss = size(pp);
pp = [pp zeros(ss(1),1)];
Eplot = E_grid;
muplot = acosd(mu_pitch_grid(1:end));

h = polarPcolor(Eplot,muplot,pp,'Ncircles',10,'Nspokes',10,'colBar',0);
cb = colorbar;
cb.Title.String = "log_{10}Ie (#e/m2/s)";
cb.Title.FontSize = 14;
cb.Position = [0.12 0.1291 0.0294 0.7775];

% h = pcolor(Eplot,(muplot),pp.');
% set(h,'EdgeColor','none');
% set(gca,'XScale','log');
% xlabel('E (eV)')
% ylabel('\theta')
% yticks([acosd(mu_pitch_grid)])
% ylim([0 90])
% cb = colorbar;
% cb.Title.String = "log_{10}Ie (#e/m2/s/eV/ster)";

% caxis([-5 2])
caxis([2 6])
%% Write Ie_total coming at top of the ionosphere in a file to be used by AURORA
clear Ie_total;
for i_t = 1:size(Ie_1_out,1)
  for i_mu = 1:9   % downward fluxes
    Ie_total{i_mu}(:,i_t) = Ie_1_out(i_t,:,i_mu);
  end
  for i_mu = 10:18  % upward fluxes
    Ie_total{i_mu}(:,i_t) = zeros(size(E_grid));
  end
end

save('Ie_incoming.mat','-v7.3','Ie_total')








%% ----------------------------------------------------------------------








%% Extract Ie_top_raw
load MIC-18streams-0.35s-3/Ie_top.mat Ie_top_raw
% E_grid = E;
Ie_out = Ie_top_raw(10:end,2:end,:); % from time t = 0.001 to 0.35,
                                    % from angle 90° to 0° (9 streams)
% f_out = zeros(size(Ie_out,2),200,100);
f_out = [];

%% Convert in the other way
% clear Ie_out DIFF
% Ie_out(:,1,:) = Ie_1(:,[end:-1:1]);
% Ie_out = permute(Ie_out,[3 2 1]);

load ./Bfield.mat
B = B(end);

% Extrapolate Ie over a finer (E,mu_pitch) grid;
HMR_E = 1;   % HOW MUCH DO YOU WANT TO REFINE E
HMR_MU = 20;   % HOW MUCH DO YOU WANT TO REFINE MU

% Refine E-grid
F = griddedInterpolant(1:701,E_grid);
E_grid_finer = F(1:(1/HMR_E):701);
dE_finer = (1/HMR_E) .* dE;
dE_finer = repelem(dE_finer,HMR_E); % or dE_finer = diff(E_grid_finer);dE_finer = [dE_finer,dE_finer(end)];
E_middle_bin_finer = E_grid_finer;

% Refine mu_pitch-grid
theta_lims2do = 90:-10:0;
mu_pitch_grid = cosd(theta_lims2do);
BeamW = [];
for iMu = (numel(mu_pitch_grid)-1):-1:1
  BeamW(iMu) = 2*pi * abs(integral(@(pa) sin(pa),acos(mu_pitch_grid(iMu)),acos(mu_pitch_grid(iMu+1))));
end
G = griddedInterpolant(1:length(theta_lims2do),theta_lims2do);  
theta_lims2do_finer = G(1:(1/HMR_MU):length(theta_lims2do));
mu_grid_finer = cos(theta_lims2do_finer*pi/180);
BeamW_finer = [];
for iMu = (numel(mu_grid_finer)-1):-1:1
  BeamW_finer(iMu) = 2*pi * abs(integral(@(pa) sin(pa),acos(mu_grid_finer(iMu)),acos(mu_grid_finer(iMu+1))));
end
mu_pitch_middle_bin_finer = mu_avg(mu_grid_finer);

vz_middle_bin = particle(index_specie).vz;

%Loop over timesteps
for i_t = 1:size(Ie_out,2)
  % Refine I
  Ie_out_finer = Ie_out(:,i_t,:)./ (BeamW([end:-1:1]).'./ sum(BeamW));
  Ie_out_finer = repelem(Ie_out_finer./(HMR_E),HMR_MU,1,HMR_E);
  Ie_out_finer = Ie_out_finer .* (BeamW_finer.' ./ sum(BeamW_finer));

  f = zeros(numel(vz_middle_bin),numel(mu_mag_middle_bin));

  v = zeros(numel(vz_middle_bin),numel(mu_mag_middle_bin));
  for ii = 1:length(mu_pitch_middle_bin_finer)
    for jj = 1:length(E_middle_bin_finer)
      % get coordinate point in (vz,mu_mag)-grid
      vz = - sqrt(2 * E_middle_bin_finer(jj) * 1.6e-19 / m) * mu_pitch_middle_bin_finer(ii);
      mu_mag = E_middle_bin_finer(jj) * 1.6e-19 ./ B * (1 - mu_pitch_middle_bin_finer(ii).^2);
      [~,index_vz] = min(abs(vz_middle_bin - vz));
      [~,index_mu_mag] = min(abs(mu_mag_middle_bin - mu_mag));
      v(index_vz,index_mu_mag) = sqrt(vz_middle_bin(index_vz).^2 + 2*B/m .* mu_mag_middle_bin(index_mu_mag));
      % convert flux and compute distribution function
      f(index_vz,index_mu_mag) = f(index_vz,index_mu_mag) + 1 ./ ...
                                v(index_vz,index_mu_mag) * ...  
                                Ie_out_finer(ii,1,jj) ./ ...
                                (dvz * dmu_mag(index_mu_mag));
    end
  end

%  Ie_total = sum((v .* f * dmu_mag.') * dvz)
% f_out(i_t,:,:) = f;
f_out = [f_out;f];
disp(['Conversion ',num2str(i_t),'/',num2str(size(Ie_out,2)),' done!'])
end

% DIFF =  abs(f - fzvzmu_in);
% Ie_total = sum((v .* DIFF * dmu_mag.') * dvz);
% Ie_total./Ie_total_1
%% Write f_out coming out from the ionosphere in a file to be used by KETCHUP

fileID = fopen('f_response.bin','w');
fwrite(fileID',f_out,'double');
fclose(fileID);
%%
% fileID = fopen('f_response.txt','w');
% fprintf(fileID','%d',f_out);
% fclose(fileID);

% writematrix(f_out,'f_response.txt','delimiter','space')

%%
% fileID = fopen('f_response.bin');
% M = fread(fileID,[30000 100],'double');
% fclose(fileID);
% size(M)

%% Plot f(vz,mu_mag)
% DIFF = log10(abs(fzvzmu_in .* (log10(fzvzmu_in) > -10) - f .* (log10(f) > -10)) ./ (fzvzmu_in .* (log10(fzvzmu_in) > -10)));
% DIFF = - (log10(fzvzmu_in) .* (log10(fzvzmu_in) > -10) - log10(f).* (log10(f) > -10));
load outp/fzvzmu0016000.mat

ii = 1;
zz = 600;

figure
c1 = jet(64);c2 = jet(256);c3 = jet(1024);
c = [c1(1:8,:);c2(32:64,:);c3(257:1024,:)];
colormap(c)

pp = log10(fzvzmustruct((ii)).f(:,:,zz));
% pp=log10(M(801:1000,:));
% pp=(DIFF(:,:));

ss = size(pp);
pp = [[pp zeros(ss(1),1)] ; zeros(1,ss(2)+1)];
vzplot = particle(ii).vzcorn;
muplot = particle(ii).mucorn;

% h = pcolor(muplot,vzplot(101:end),pp);
h = pcolor(muplot,vzplot,pp);

axis([min(particle(ii).mucorn) ...
      max(particle(ii).mucorn) ...
      min(particle(ii).vzcorn) ...
      max(particle(ii).vzcorn)])
set(h,'EdgeColor','none');
cb = colorbar;
cb.Title.String = "log_{10}(f(v_z,\mu_{mag})";
cax = caxis;
caxis([-10 cax(2)]) 
xlabel('\mu_{mag}')
ylabel('v_z')
title(['z = ',num2str(z(zz),3)])

% xlim([0 2e-13])   %ionospheric
% ylim([-2e7 2e7])  %ionospheric

xlim([0 9e-11]) %magnetospheric
ylim([-1e8 1e8])   %magnetospheric
caxis([-10 2])
%% Plot f(vz,mu_mag) at the ionospheric boundary (BC !)
% Need that the fR0000000s03.ketchup.dat file has been copied in the
% folder, and that the extract_BC.m script has been used to convert it
% (the .dat file) into a .mat file.
load fBC_right_s03.mat
ii = 3;
 
figure
c1 = jet(64);c2 = jet(256);c3 = jet(1024);
c = [c1(1:8,:);c2(32:64,:);c3(257:1024,:)];
colormap(c)
pp = log10(fBC(ii).f);

ss = size(pp);
pp = [[pp zeros(ss(1),1)] ; zeros(1,ss(2)+1)];
vzplot = particle(ii).vzcorn;
muplot = particle(ii).mucorn;
h = pcolor(muplot,vzplot,pp);
set(h,'EdgeColor','none');
cb = colorbar;
cb.Title.String = "log_{10}(f(v_z,\mu_{mag})";
cax = caxis;
caxis([-10 cax(2)]) 
xlabel('\mu_{mag}')
ylabel('v_z')
title(['z = ionospheric BC'])

xlim([0 2e-13])   %ionospheric
ylim([-2e7 2e7])  %ionospheric
%% Contour plot of pitch angles
mu_pitch = [];
for ii = 1:length(vz_grid)
  for jj = 1:length(mu_mag_grid)
        mu_pitch(ii,jj) = cos(atan(sqrt(2*B/m * mu_mag_grid(jj)) / vz_grid (ii)));
  end
end
[X,Y] = meshgrid(mu_mag_grid,vz_grid);
figure(5)
contour(X,Y,acosd(mu_pitch),180-theta_lims2do-5,'ShowText','on','linewidth',1)
xlim([0 6e-11])
ylim([0 1e8])

