% This script converts the distribution function extracted from the
% MI-model and expressed as a function of vz and mu_mag into a distribution
% function expressed as a function of E and mu_pitch that can be used as
% input in a ionosphere response model.

%% Import data
inputb6;
load fzvzmu3000000.mat
load Bfield.mat
%% Calculate flux (/eV/ster) for specie 1 (magnetospheric e-)
index_specie = 1;
m = particle(index_specie).mass;
zz = 531;      % index for the altitude at top of ionosphere
B = B(end);   % magnetic field at top of ionosphere

% Extraction of the (vz,mu_mag)-grid
vz_middle_bin = particle(index_specie).vz(particle(index_specie).vz >= 0);                 % taking only downward going e-
vz_grid = particle(index_specie).vzcorn(particle(index_specie).vzcorn >= 0);    % taking only downward going e-
dvz = particle(index_specie).dvz;
mu_mag_middle_bin = particle(index_specie).mu;
mu_mag_grid = particle(index_specie).mucorn;
dmu_mag = particle(index_specie).dmu;
% Extraction of the distr funct f(v,mu_mag)
fzvzmu_in = fzvzmustruct(index_specie).f(particle(index_specie).vz >= 0,:,zz);  % taking only downward going e-

% Extrapolate f over a finer (vz,mu_mag) grid;
HMR_VZ = 50;   % HOW MUCH DO YOU WANT TO REFINE VZ
HMR_MU = 50;   % HOW MUCH DO YOU WANT TO REFINE MU
% Refine vz-grid
F = griddedInterpolant(1:251,vz_grid);
vz_grid_finer = F(1:(1/HMR_VZ):251);
dvz_finer = (1/HMR_VZ) * dvz ;
vz_middle_bin_finer = vz_grid_finer(2:end) - 0.5 * dvz_finer;
% Refine mu_mag-grid
G = griddedInterpolant(1:51,mu_mag_grid);
mu_grid_finer = G(1:(1/HMR_MU):51);
dmu_mag_finer = (1/HMR_MU) .* dmu_mag;
dmu_mag_finer = repelem(dmu_mag_finer,HMR_MU);
mu_mag_middle_bin_finer = mu_grid_finer(2:end) - 0.5 .* dmu_mag_finer;
% Refine f
fzvzmu_in_finer = repelem(fzvzmu_in,HMR_VZ,HMR_MU);

% Make E-grid
dEfcn = @(X,DEinitial,DEfinal,C,X0) DEinitial+(1+tanh(C*(X-X0)))/2*DEfinal;     % from setup4etrpd10streams.m
E_grid = cumsum(dEfcn(0:2000,0.15,11.5,0.05,80))+1.9;
E_grid = E_grid;
dE = diff(E_grid);
E_middle_bin = E_grid(1:(end-1)) + 0.5 * dE;
% Make mu_pitch-grid
% dtheta = [0 10 20 30 20 10 10 20 30 20 10];
% theta_lims2do = [180 170 150 120 100 90];
theta_lims2do = [180 170 150 130 125 120 115 110 105 100 95 90];
[~,~,BeamW,mu_pitch_grid] = e_scattering_result_finder(theta_lims2do,AURORA_root_directory);
mu_pitch_grid = abs(mu_pitch_grid);  % taking downward going e-
BeamW = 2*pi*BeamW;                  % taking downward going e-
mu_pitch_middle_bin = mu_avg(mu_pitch_grid);

% mu_pitch_grid = cos((0:10:90)*pi/180);
% dmu_pitch = abs(diff(mu_pitch_grid));
% mu_pitch_middle_bins = mu_avg(mu_pitch_grid);


% RELATIVISTIC PARTICLES
% The magnetospheric electrons are considered as relativistic particles in
% the ketchup code. This means that the magnetic moment mu_mag' associated
% to the particles has been corrected with the lorentz factor gamma :
%               mu_mag' = gamma * mu_mag 
% For our conversion here to (E,mu_pitch) grid, we need to convert mu_mag'
% back into mu_mag
% c = 3e8;
% gamma = @(v,mu) sqrt( (1 + 2*B*mu/c^2/m) ./ (1 - v^2/c^2) );
gamma = @(v,mu) 1;


% vz_middle_bin = vz_middle_bin_finer;
% dvz = dvz_finer;
% mu_mag_middle_bin = mu_mag_middle_bin_finer;
% dmu_mag = dmu_mag_finer;
% fzvzmu_in = fzvzmu_in_finer;
Ie_1 = zeros(numel(E_middle_bin),numel(mu_pitch_middle_bin));

for ii = 1:length(vz_middle_bin_finer)
  for jj = 1:length(mu_mag_middle_bin_finer)
    % get coordinate point in (E,mu_pitch)-grid
    E = 0.5 * m * (vz_middle_bin_finer(ii).^2 + 2*B/m * 1/gamma(vz_middle_bin_finer(ii),mu_mag_middle_bin_finer(jj)) .* mu_mag_middle_bin_finer(jj)) ./ 1.6e-19;
    mu_pitch = cos(atan(sqrt(2*B/m * 1/gamma(vz_middle_bin_finer(ii),mu_mag_middle_bin_finer(jj)) * mu_mag_middle_bin_finer(jj)) / vz_middle_bin_finer (ii)));
    [~,index_pitch] = min(abs(mu_pitch_middle_bin - mu_pitch));
    [~,index_energy] = min(abs(E_middle_bin - E));
    % convert distribution function and compute flux (/eV/ster)
    Ie_1(index_energy,index_pitch) = Ie_1(index_energy,index_pitch) + ...
                                vz_middle_bin_finer(ii) * fzvzmu_in_finer(ii,jj) * ...
                                dvz_finer * dmu_mag_finer(jj) / ...
                                (dE(index_energy) * BeamW(index_pitch));
  end
end

Ie_total_new1 = (Ie_1 * BeamW.').' * dE.';                % *dmu_pitch *dE [#e/m2/s]
Ie_total_old1 = sum(((vz_middle_bin_finer.' .* fzvzmu_in_finer) * dmu_mag_finer.') * dvz_finer); % *dmu_mag * dvz [#e/m2/s]

%% Calculate flux (/eV/ster) for specie 3 (ionospheric e-)

index_specie = 3;
m = particle(index_specie).mass;
zz = Nz;      % index for the altitude at top of ionosphere
B = B(end);   % magnetic field at top of ionosphere

% Extraction of the (vz,mu_mag)-grid
vz_middle_bin = particle(index_specie).vz(particle(index_specie).vz >= 0);                 % taking only downward going e-
vz_grid = particle(index_specie).vzcorn(particle(index_specie).vzcorn >= 0);    % taking only downward going e-
dvz = particle(index_specie).dvz;
mu_mag_middle_bin = particle(index_specie).mu;
mu_mag_grid = particle(index_specie).mucorn;
dmu_mag = particle(index_specie).dmu;
% Extraction of the distr funct f(v,mu_mag)
fzvzmu_in = fzvzmustruct(index_specie).f(particle(index_specie).vz >= 0,:,zz);  % taking only downward going e-

% Extrapolate f over a finer (vz,mu_mag) grid;
HMR_VZ = 10;   % HOW MUCH DO YOU WANT TO REFINE VZ
HMR_MU = 10;   % HOW MUCH DO YOU WANT TO REFINE MU
% Refine vz-grid
F = griddedInterpolant(1:251,vz_grid);
vz_grid_finer = F(1:(1/HMR_VZ):251);
dvz_finer = (1/HMR_VZ) * dvz ;
vz_middle_bin_finer = vz_grid_finer(2:end) - 0.5 * dvz_finer;
% Refine mu_mag-grid
G = griddedInterpolant(1:51,mu_mag_grid);
mu_grid_finer = G(1:(1/HMR_MU):51);
dmu_mag_finer = (1/HMR_MU) .* dmu_mag;
dmu_mag_finer = repelem(dmu_mag_finer,HMR_MU);
mu_mag_middle_bin_finer = mu_grid_finer(2:end) - 0.5 .* dmu_mag_finer;
% Refine f
fzvzmu_in_finer = repelem(fzvzmu_in,HMR_VZ,HMR_MU);

% Make E-grid
dEfcn = @(X,DEinitial,DEfinal,C,X0) DEinitial+(1+tanh(C*(X-X0)))/2*DEfinal;     % from setup4etrpd10streams.m
E_grid = cumsum(dEfcn(0:937,0.15,11.5,0.05,80))+1.9;
E_grid = E_grid;
dE = diff(E_grid);
E_middle_bin = E_grid(1:(end-1)) + 0.5 * dE;
% Make mu_pitch-grid

dtheta = [0 10 20 30 20 10 10 20 30 20 10];
theta_lims2do = 180-cumsum(dtheta);
[~,~,BeamW,mu_pitch_grid] = e_scattering_result_finder(theta_lims2do,AURORA_root_directory);
mu_pitch_grid = abs(mu_pitch_grid(1:6));  % taking downward going e-
BeamW = BeamW(1:5);                       % taking downward going e-
mu_pitch_middle_bin = mu_avg(mu_pitch_grid);
% mu_pitch_grid = cos((0:10:90)*pi/180);
% mu_pitch_grid = cos(theta_lims2do);
% dmu_pitch = abs(diff(mu_pitch_grid));



% RELATIVISTIC PARTICLES
% The magnetospheric electrons are considered as relativistic particles in
% the ketchup code. This means that the magnetic moment mu_mag' associated
% to the particles has been corrected with the lorentz factor gamma :
%               mu_mag' = gamma * mu_mag 
% For our conversion here to (E,mu_pitch) grid, we need to convert mu_mag'
% back into mu_mag
c = 3e8;
gamma = @(v,mu) sqrt( (1 + 2*B*mu/c^2/m) ./ (1 - v^2/c^2) );


vz_middle_bin = vz_middle_bin_finer;
dvz = dvz_finer;
mu_mag_middle_bin = mu_mag_middle_bin_finer;
dmu_mag = dmu_mag_finer;
fzvzmu_in = fzvzmu_in_finer;
Ie_3 = zeros(numel(E_middle_bin),numel(mu_pitch_middle_bin));

for ii = 1:length(vz_middle_bin)
  for jj = 1:length(mu_mag_middle_bin)
    % get coordinate point in (E,mu_pitch)-grid
    E = 0.5 * m * (vz_middle_bin(ii).^2 + 2*B/m * 1/gamma(vz_middle_bin(ii),mu_mag_middle_bin(jj)) .* mu_mag_middle_bin(jj)) ./ 1.6e-19;
    mu_pitch = cos(atan(sqrt(2*B/m * 1/gamma(vz_middle_bin(ii),mu_mag_middle_bin(jj)) * mu_mag_middle_bin(jj)) / vz_middle_bin (ii)));
    [~,index_pitch] = min(abs(mu_pitch_middle_bin - mu_pitch));
    [~,index_energy] = min(abs(E_middle_bin - E));
    % convert distribution function and compute flux (/eV/ster)
    Ie_3(index_energy,index_pitch) = Ie_3(index_energy,index_pitch) + ...
                                vz_middle_bin(ii) * fzvzmu_in(ii,jj) * ...
                                dvz * dmu_mag(jj) / ...
                                (dE(index_energy) * BeamW(index_pitch));
  end
end

Ie_total_new3 = (Ie_3 * dmu_pitch.').' * dE.';                % *dmu_pitch *dE [#e/m2/s]
Ie_total_old3 = sum(((vz_middle_bin.' .* fzvzmu_in) * dmu_mag.') * dvz); % *dmu_mag * dvz [#e/m2/s]
%% plot Ie
Ie_plot = Ie_1;
% Ie_plot = Ie_3;
% Ie_plot = Ie_1 + Ie_3; 

figure()
c1=jet(64);c2=jet(256);c3=jet(1024);
c=[c1(1:8,:);c2(32:64,:);c3(257:1024,:)];
colormap(c)
pp=log10(Ie_plot);
ss=size(pp);
pp=[[pp zeros(ss(1),1)];zeros(1,ss(2)+1)];
Eplot=E_grid;
muplot=acosd(mu_pitch_grid(1:end));
h = pcolor(Eplot,muplot,pp.');
cb = colorbar;
yticks([acosd(mu_pitch_grid)])
ylim([0 90])
% view(2);grid on;shading faceted %shading flat

% axis([min(E_grid) ...
%       max(E_grid) ...
%       min(muplot) ...
%       max(muplot)])
% fig = get(gca)
set(h,'EdgeColor','none');
set(gca,'XScale','log');
xlabel('E (eV)')
ylabel('\theta')
cb.Title.String = "log_{10}Ie (#e/m2/s/eV/ster)";
caxis([-10 1])
% Eplot=E_middle_bins;
% muplot=acosd(mu_pitch_grid(2:end))-5;
% imagesc(Eplot,muplot,log10(Ie).')
% colorbar
% yticks([acosd(mu_pitch_grid)])
% ylim([0 90])
%% Compute total flux
Ie_total_new = Ie_total_new1 %+ Ie_total_new3  % [#e/m2/s]
Ie_total_old = Ie_total_old1 %+ Ie_total_old3  % [#e/m2/s
%% Convert in the other way
% Extrapolate Ie over a finer (E,mu_pitch) grid;
HMR_E = 1;   % HOW MUCH DO YOU WANT TO REFINE E
HMR_MU = 50;   % HOW MUCH DO YOU WANT TO REFINE MU
% Refine E-grid
F = griddedInterpolant(1:2001,E_grid);
E_grid_finer = F(1:(1/HMR_E):2001);
dE_finer = (1/HMR_E) .* dE;
dE_finer = repelem(dE_finer,HMR_E);
E_middle_bin_finer = E_grid_finer(2:end) - 0.5 .* dE_finer;
% Refine mu_pitch-grid
G = griddedInterpolant(1:length(theta_lims2do),theta_lims2do);  
theta_lims2do_finer = G(1:(1/HMR_MU):length(theta_lims2do));
[~,~,BeamW_finer,mu_grid_finer] = e_scattering_result_finder(theta_lims2do_finer,AURORA_root_directory);
mu_grid_finer = abs(mu_grid_finer);
BeamW_finer = 2*pi*BeamW_finer;
mu_pitch_middle_bin_finer = mu_avg(mu_grid_finer);
% Refine I
Ie_1_finer = repelem(Ie_1,HMR_E,HMR_MU);

f = zeros(numel(vz_middle_bin),numel(mu_mag_middle_bin));
for ii = 1:length(E_middle_bin_finer)
  for jj = 1:length(mu_pitch_middle_bin_finer)
    % get coordinate point in (vz,mu_mag)-grid
    vz = sqrt(2 * E_middle_bin_finer(ii) * 1.6e-19 / m) * mu_pitch_middle_bin_finer(jj);
    mu_mag = E_middle_bin_finer(ii) * 1.6e-19 ./ B * (1 - mu_pitch_middle_bin_finer(jj).^2);
    [~,index_vz] = min(abs(vz_middle_bin - vz));
    [~,index_mu_mag] = min(abs(mu_mag_middle_bin - mu_mag));
    % convert flux and compute distribution function
    f(index_vz,index_mu_mag) = f(index_vz,index_mu_mag) + 1 ./ vz_middle_bin(index_vz) * ...
                               Ie_1_finer(ii,jj) * dE_finer(ii) * BeamW_finer(jj) ./ ...
                               (dvz * dmu_mag(index_mu_mag));
  end
end

Ie_total = sum(((vz_middle_bin.' .* f) * dmu_mag.') * dvz)

%% Plot f(vz,mu_mag)
finitemass=[];
for ii=1:Nspecies
    if ~isnan(particle(ii).mass) & ~isinf(particle(ii).mass)
        finitemass=[finitemass ii];
    end
end
ii = 1;
zz = 531;

figure
c1=jet(64);c2=jet(256);c3=jet(1024);
c=[c1(1:8,:);c2(32:64,:);c3(257:1024,:)];
colormap(c)
pp=log10(fzvzmustruct((ii)).f(:,:,zz));
% pp=log10(f(:,:));
ss=size(pp);
pp=[[pp zeros(ss(1),1)];zeros(1,ss(2)+1)];
vzplot=particle(finitemass(ii)).vzcorn + ...
       fzvzmustruct(finitemass(ii)).ivzoffset(zz) * ...
       particle(finitemass(ii)).dvz;
muplot=particle(finitemass(ii)).mucorn;
% h = surf(muplot,vzplot(251:end),pp);
h = surf(muplot,vzplot(251:end),pp(251:end,:));
view(2);grid on;shading faceted %shading flat
axis([min(particle(finitemass(ii)).mucorn) ...
      max(particle(finitemass(ii)).mucorn) ...
      (min(particle(finitemass(ii)).vzcorn) + ...
       fzvzmustruct(finitemass(ii)).ivzoffset(zz) * ...
       particle(finitemass(ii)).dvz) ...
      (max(particle(finitemass(ii)).vzcorn) + ...
       fzvzmustruct(finitemass(ii)).ivzoffset(zz) * ...
       particle(finitemass(ii)).dvz)])
% set(gca,'fontname','utopia','fontsize',14)
set(h,'EdgeColor','none');
cb = colorbar;
cax=caxis;
caxis([-10 cax(2)]) 
xlim([0 6e-11])
xlabel('\mu_{mag}')
ylabel('v_z')
cb.Title.String = "log_{10}(f(v_z,\mu_{mag})";
