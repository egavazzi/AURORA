function [Pmu2mup,theta2beamW,BeamW,mu_lims] = e_scattering_result_lister(AURORA_root_directory)
% E_SCATTERING_RESULT_FINDER - scattering-matrix-loading/calculating
%   E_SCATTERING_RESULT_FINDER tries to load the beam-2-beam
%   scattering matrices from the E_scatterings directory matching
%   the recuested pitch-angle limits (theta_lims2do) in the
%   AURORA_roo_directory. If no file with matching
%   pitch-angle-limits are found new calculations are done with
%   e_scattering_beamdistribution, and are saved in the
%   E_scattering directory. This way the unreasonable slow
%   calculations of these scattering-matrices will only be done
%   once for each pitch-angle-stream configuration and subsequent
%   runs will utilize the correct pre-calculated matrices.
% 
% Calling:
%   [Pmu2mup,theta2beamW,BeamW,mu_lims] = e_scattering_result_finder(theta_lims2do,AURORA_root_directory)
% Input:
%  theta_lims2do - pitch-angle-limit array (deg) [1 x (n_mu+1)]
%                  double array
%  AURORA_root_directory - full path to the AURORA-root directory,
%                  automatically set with the add_AURORA script
% Output:
%  Pmu2mup     - double array [n_dirs x n_dirs x n_beams] with the
%                probability distribution for an electron ending up
%                in a beam starting at pitch-angle theta0 and
%                scattering an angle theta1.
%  theta2beamW - Relative weighting matrix, double matrix [n_beams
%                x n_dirs], with the relative contribution from
%                within each beam, assuming isotropic pich-angle
%                distribution within each beam.
%  BeamW       - pitch-angle-stream solid angle array (ster) double
%                array [1 x n_mu]
%  mu_lims     - cos(theta_lims)
% 
% SEE ALSO: e_scattering_beamdistribution

%  Copyright © Bjorn Gustavsson 20191122, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later

% All pitch-angle-scattering files with possible matching theta_lims:
e_s_b_files = dir(fullfile(AURORA_root_directory,'E_scatterings'));

char(e_s_b_files(:).name)
for i1 = 1:numel(e_s_b_files),
  if ~ e_s_b_files(i1).isdir
    %try
      load(fullfile(AURORA_root_directory,...
                    'E_scatterings',...
                    e_s_b_files(i1).name),...
           'theta_lims')
      % then we have found a match, so load
      theta_str = sprintf('%3.0f, ',theta_lims);
      fprintf('pitch-angle limits: %s (deg) ',theta_str(1:end-2))
      fprintf(' file: %s\n',e_s_b_files(i1).name)
    %catch
      
      % do nothing, bad practise perhaps, but I cannot be bovvered
      % to write test-exceptionings for this task of attempting to
      % load mat-files in a directory and check if the
      % possible/probable variable theta_lims matches theta_lims2do.
      % It is simple enough to just use try and see if we get a
      % match, then if there is none we still know that some work
      % needs to be done below...
    %end
  end
end
