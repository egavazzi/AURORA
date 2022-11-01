function Ie_total = conversion_M_to_I(path_to_vlasov_initial_file,path_to_vlasov_simulation_results,index_specie,HMR_VZ,HMR_MU,E_max,theta_lims2do,first_run)


cd(path_to_vlasov_simulation_results)
%% Import data
inputb6; 
load outp/Bfield.mat;


%% Calculate flux ([#e/m2/s]) for specie 1 (magnetospheric e-)
% Extracting fzvzmu_in
if first_run
  load(path_to_vlasov_initial_file)
  f_in(1,:,:) = fzvzmustruct(1).f(particle(1).vz >= 0,:,end-1);
  disp(['Loading fzvzmu_in (1/1)'])
else
  % displays important informations to the user about the total time of
  % simulation and time intervals between data points
  disp(['Simulation time : ',num2str(dt*Niter),'s'])
  disp(['Time between distr. data points :',num2str(dt*dump_period_distr),'s'])
  disp(['--------------------'])
  cd ./f_vlasov_raw/
  dDir = dir;
  f_in = zeros(numel(dDir)-1,100,100);
  load(path_to_vlasov_initial_file)
  f_in(1,:,:) = fzvzmustruct(1).f(particle(1).vz >= 0,:,end-1);
  step = 1; steps = size(f_in,1);
  disp(['Loading fzvzmu_in (',num2str(step),'/',num2str(steps),')'])
  for i1 = 3:numel(dDir)
    step = i1-1;
    load(dDir(i1).name);
    f_in(i1-1,:,:) = fzvzmustruct(1).f(particle(1).vz >= 0,:,end-1);    % taking only downward going e-
    disp(['Loading fzvzmu_in (',num2str(step),'/',num2str(steps),')'])
  end
cd ../
end

m = particle(index_specie).mass;   
B = B(end-1); % magnetic field at top of ionosphere
% Extraction of the (vz,mu_mag)-grid
vz_middle_bin = particle(index_specie).vz(particle(index_specie).vz >= 0);      % taking only downward going e-
vz_grid = particle(index_specie).vzcorn(particle(index_specie).vzcorn >= 0);    % taking only downward going e-
dvz = particle(index_specie).dvz;
mu_mag_middle_bin = particle(index_specie).mu;
mu_mag_grid = particle(index_specie).mucorn;
dmu_mag = particle(index_specie).dmu;

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
dEfcn = @(X,DEinitial,DEfinal,C,X0) DEinitial+(1+tanh(C*(X-X0)))/2*DEfinal;
E_grid = cumsum(dEfcn(0:2000,0.15,11.5,0.05,80))+1.9;
[~,iE_max] = min(abs(E_grid-E_max)); 
E_middle_bin = E_grid(1:iE_max);
dE = diff(E_grid);dE = [dE,dE(end)];
% Make mu_pitch-grid
mu_pitch_grid = cosd(theta_lims2do);
BeamW = [];
for iMu = (numel(mu_pitch_grid)-1):-1:1
  BeamW(iMu) = 2*pi * abs(integral(@(pa) sin(pa),acos(mu_pitch_grid(iMu)),acos(mu_pitch_grid(iMu+1))));
end
mu_pitch_middle_bin = mu_avg(mu_pitch_grid);

Ie_in = zeros(size(f_in,1),length(E_middle_bin),length(mu_pitch_middle_bin));
step = 1; steps = size(f_in,1);
for i_t = 1:steps
  step = i_t;
  disp(['Converting f into I (',num2str(step),'/',num2str(steps),')'])
  % Extraction of the distr funct f(v,mu_mag)
  fzvzmu_in = f_in(i_t,:,:);
  fzvzmu_in = squeeze(fzvzmu_in);
  % Extrapolate f over the finer (vz,mu_mag) grid;
  fzvzmu_in_finer = repelem(fzvzmu_in,HMR_VZ,HMR_MU);
  
  Ie_temp = zeros(numel(E_middle_bin),numel(mu_pitch_middle_bin));
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
      Ie_temp(index_energy,index_pitch) = Ie_temp(index_energy,index_pitch) + ...
                                  v(ii,jj) * fzvzmu_in_finer(ii,jj) * ...
                                  dvz_finer * dmu_mag_finer(jj);
    end
  end
  Ie_in(i_t,:,:) = Ie_temp;
end

clear Ie_total;
for i_t = 1:size(Ie_in,1)
  % downward fluxes
  for i_mu = 1:(length(theta_lims2do)-1)
    Ie_total{i_mu}(:,i_t) = Ie_in(i_t,:,i_mu);
  end
  % upward fluxes
  for i_mu = (length(theta_lims2do)-1):2*(length(theta_lims2do)-1)
    Ie_total{i_mu}(:,i_t) = zeros(size(E_middle_bin));
  end
end