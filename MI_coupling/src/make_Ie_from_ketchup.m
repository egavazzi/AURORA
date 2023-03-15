function [Ie_DL,z_DL] = make_Ie_from_ketchup(path_to_vlasov_initial_file, ...
                                              path_to_vlasov_simulation_results, ...
                                              index_specie, zmax_to_extract,h_atm_iono,HMR_VZ,HMR_MU,iE_max,mu_lims)
% This script extracts and converts the function distributions from ketchup into fluxes of
% electrons Ie.
% It looks through a directory with all the fzvzmu***.mat files from a simulation. Each file 
% should correspond to one time step and contain the function distribution over all the space 
% points of the simulation.

%% Import data
run(fullfile(path_to_vlasov_simulation_results,'inputb6.m'));
load(fullfile(path_to_vlasov_simulation_results,'outp/Bfield.mat'))
Nvz = particle(index_specie).Nvz;
Nmu = particle(index_specie).Nmu;

% Extracting fzvzmu_in
dDir = dir(fullfile(path_to_vlasov_simulation_results,'f_vlasov_raw'));
f_in_all = zeros(numel(dDir)-1,Nvz,Nmu,Nz);
load(path_to_vlasov_initial_file)
f_in_all(1,:,:,:) = fzvzmustruct(1).f(:,:,:);
step = 1; steps = size(f_in_all,1);
disp(['Loading fzvzmu (',num2str(step),'/',num2str(steps),')'])
for i1 = 3:numel(dDir)
  step = i1-1;
  load(fullfile(path_to_vlasov_simulation_results,'f_vlasov_raw/',dDir(i1).name));
  f_in_all(i1-1,:,:,:) = fzvzmustruct(1).f(:,:,:);
  disp(['Loading fzvzmu (',num2str(step),'/',num2str(steps),')'])
end


m = particle(index_specie).mass;

% Convert z into km ABOVE the ionosphere 
% ignoring the two points in ketchup that overlap with the ionospheric simulation)
% z_DL = z(end) - z(end-1:-1:1) + h_atm_iono(end);
% NOT ignoring the two points in ketchup that overlap with the ionospheric simulation)
z_DL = z(end) - z(end:-1:1) + h_atm_iono(end-1); 
z_DL = z_DL/1e3;

% and find the index for max altitude to extract data from
[~, i_zmax] = min(abs(z_DL - zmax_to_extract));
z_DL = z_DL(1:i_zmax);

%% Loop over the altitudes to extract ([#e/m2/s])
step = 1; steps = i_zmax;
Ie_DL = zeros(i_zmax*(numel(mu_lims)-1), size(f_in_all,1), iE_max);
for i_z = 1:steps
  
  % Convert index to match matrices from ketchup. 
%   zz = Nz-i_z-1; % ignoring the two points overlaping with the ionospheric simulation.
  zz = Nz-i_z+1; % NOT ignoring the two points overlaping with the ionospheric simulation.

  step = i_z;
  disp(['Converting f into I over altitudes (',num2str(step),'/',num2str(steps),')'])
    
  B_z = B(zz);   % magnetic field at altitude zz
  f_in = squeeze(f_in_all(:,:,:,zz)); % distr. func. at altitude zz

  % Extraction of the (vz,mu_mag)-grid
  vz_middle_bin = particle(index_specie).vz;      % taking only downward going e-
  vz_grid = particle(index_specie).vzcorn;        % taking only downward going e-
  dvz = particle(index_specie).dvz;
  mu_mag_middle_bin = particle(index_specie).mu;
  mu_mag_grid = particle(index_specie).mucorn;
  dmu_mag = particle(index_specie).dmu;

  % Refine vz-grid
  F = griddedInterpolant(1:201,vz_grid);
  vz_grid_finer = F(1:(1/HMR_VZ):201);
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
  E_grid = cumsum(dEfcn(0:2000,0.15,11.5,0.05,80))+1.9;
  E_middle_bin = E_grid(1:iE_max);
  dE = diff(E_grid);dE = [dE,dE(end)];
  % Make mu_pitch-grid
  mu_pitch_grid = mu_lims;
  BeamW = [];
  for iMu = (numel(mu_pitch_grid)-1):-1:1
    BeamW(iMu) = 2*pi * abs(integral(@(pa) sin(pa),acos(mu_pitch_grid(iMu)),acos(mu_pitch_grid(iMu+1))));
  end
  mu_pitch_middle_bin = mu_avg(mu_pitch_grid);

  for i_t = 1:size(f_in,1);
    % Extraction of the distr funct f(v,mu_mag)
    fzvzmu_in = f_in(i_t,:,:);
    fzvzmu_in = squeeze(fzvzmu_in);
    % Refine f
    fzvzmu_in_finer = repelem(fzvzmu_in,HMR_VZ,HMR_MU);

    Ie_temp = zeros(numel(E_middle_bin),numel(mu_pitch_middle_bin));
    v = zeros(length(vz_middle_bin),length(mu_mag_middle_bin));
    for ii = 1:length(vz_middle_bin)
      for jj = 1:length(mu_mag_middle_bin)
        v(ii,jj) = sqrt(vz_middle_bin(ii).^2 + 2*B_z/m .* mu_mag_middle_bin(jj));
      end
    end
    v = repelem(v,HMR_VZ,HMR_MU);

    for ii = 1:length(vz_middle_bin_finer)
      for jj = 1:length(mu_mag_middle_bin_finer)
        % get coordinate point in (E,mu_pitch)-grid
        E = 0.5 * m * (vz_middle_bin_finer(ii).^2 + 2*B_z/m .* mu_mag_middle_bin_finer(jj)) ./ 1.6e-19;
        mu_pitch = -sign(vz_middle_bin_finer(ii)) * cos(atan(sqrt(2*B_z/m * mu_mag_middle_bin_finer(jj)) / vz_middle_bin_finer(ii)));
        [~,index_pitch] = min(abs(mu_pitch_middle_bin - mu_pitch));
        [~,index_energy] = min(abs(E_middle_bin - E));
        % convert distribution function and compute flux
        Ie_temp(index_energy,index_pitch) = Ie_temp(index_energy,index_pitch) + ...
                                    v(ii,jj) * fzvzmu_in_finer(ii,jj) * ...
                                    dvz_finer * dmu_mag_finer(jj);
      end
    end
    
    for ii = 1:length(mu_pitch_middle_bin)
      Ie_DL(i_zmax*(ii-1) + i_z,i_t,:) = Ie_temp(:,ii);
    end
    
  end
  
end