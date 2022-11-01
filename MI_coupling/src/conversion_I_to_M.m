function [f_out,f_out2] = conversion_I_to_M(path_to_vlasov_initial_file,path_to_vlasov_simulation_results,path_to_aurora_simulation,index_specie,HMR_E,HMR_MU,E_max,theta_lims2do)

cd(path_to_vlasov_simulation_results)
%% Import data
inputb6; 
load outp/Bfield.mat;
load(path_to_vlasov_initial_file);

cd(path_to_aurora_simulation)
dDir = dir;
for i1 = 3:numel(dDir)
	if dDir(i1).isdir
		cd(dDir(i1).name)
		load Ie_MI_top.mat Ie_top_raw Ie_top2_raw E
	end
end
cd(path_to_aurora_simulation)

% Ie_out = Ie_top_raw(size(Ie_top_raw,1)/2+1:end,1:end,:); 	% from time t = 0 to 0.35, from angle 90° to 0°
% Ie_out2 = Ie_top2_raw(size(Ie_top2_raw,1)/2+1:end,1:end,:);  % from time t = 0 to 0.35, from angle 90° to 0°
Ie_out = Ie_top_raw;
Ie_out2 = Ie_top2_raw;
                                    													
% Resize E_grid																															
[~,iE_max] = min(abs(E-E_max));
E_grid = E(1:iE_max);
dE = diff(E_grid);dE = [dE,dE(end)];

% Extraction of the (vz,mu_mag)-grid
m = particle(index_specie).mass;
B = B(end-1); % magnetic field at top of ionosphere
vz_middle_bin = particle(index_specie).vz;
dvz = particle(index_specie).dvz;
mu_mag_middle_bin = particle(index_specie).mu;
dmu_mag = particle(index_specie).dmu;

% Refine E-grid
F = griddedInterpolant(1:length(E_grid),E_grid);
E_grid_finer = F(1:(1/HMR_E):length(E_grid)+1);
dE_finer = (1/HMR_E) .* dE;
dE_finer = repelem(dE_finer,HMR_E); % or dE_finer = diff(E_grid_finer);dE_finer = [dE_finer,dE_finer(end)];
E_middle_bin_finer = E_grid_finer(2:end) - 0.5*dE_finer;
% Refine mu_pitch-grid
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

f_out = [];
f_out2 = [];
step = 1; steps = size(Ie_out,2);
for i_t = 1:steps
	step = i_t;
	disp(['Converting I into f (',num2str(step),'/',num2str(steps),')'])
  % Extrapolate Ie over the finer (E,mu_pitch) grid;
  Ie_out_finer = Ie_out(:,i_t,:)./ (BeamW([end:-1:1]).'./ sum(BeamW));
  Ie_out_finer = repelem(Ie_out_finer./(HMR_E),HMR_MU,1,HMR_E);
  Ie_out_finer = Ie_out_finer .* (BeamW_finer.' ./ sum(BeamW_finer));

	Ie_out2_finer = Ie_out2(:,i_t,:)./ (BeamW([end:-1:1]).'./ sum(BeamW));
  Ie_out2_finer = repelem(Ie_out2_finer./(HMR_E),HMR_MU,1,HMR_E);
  Ie_out2_finer = Ie_out2_finer .* (BeamW_finer.' ./ sum(BeamW_finer));

  f = zeros(numel(vz_middle_bin),numel(mu_mag_middle_bin));
	f2 = zeros(numel(vz_middle_bin),numel(mu_mag_middle_bin));
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
			f2(index_vz,index_mu_mag) = f(index_vz,index_mu_mag) + 1 ./ ...
                                v(index_vz,index_mu_mag) * ...  
                                Ie_out2_finer(ii,1,jj) ./ ...
                                (dvz * dmu_mag(index_mu_mag));
    end
  end
	f_out = [f_out;f];
	f_out2 = [f_out2;f];
end