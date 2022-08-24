%% Electron-flux-at-top-of-ionosphere-extraction
% This script extracts the electron-fluxes at the top of the
% ionosphere for all pitch-angle-streams. Same as
% make_all_Ie_top_MI but extracts also z_max-1 in addition to
% z_max.

%for i2 = 1:numel(current_folder)
  cd(current_folder)
  dDir = dir;
  for i1 = 3:numel(dDir)
    if dDir(i1).isdir
      cd(dDir(i1).name)
      CD = pwd;
      try
      [t,h_atm,E,mu_lims,Ie_ZTE,mu_scatterings] = Ie_ztE_loader({'.'});
      dE = diff(E);
      dE = dE([1:end,end]);
      BeamW = mu_scatterings{3};
      nZ = numel(h_atm);

      Ie_top_raw = Ie_ZTE(nZ:nZ:end,:,:);
      Ie_top2_raw = Ie_ZTE(nZ-1:nZ:end,:,:);
  
      save('Ie_MI_top.mat','Ie_top_raw','Ie_top2_raw','E','dE','BeamW','t')
      clear Ie_top_raw Ie_top2_raw
      fprintf(':::Processed OK: %s\n',CD)
      catch
        fprintf('Everythings not right in directory: %s\n',CD)
      end
    end
  end
%end

cd(current_folder)