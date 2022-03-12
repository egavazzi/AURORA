function [S2ndO,S2ndO2,S2ndN2] = e_cascading_result_finder(E_2do,AURORA_root_directory)
% E_CASCADING_RESULT_FINDER - cascading-matrix-loading/calculating
%   E_CASCADING_RESULT_FINDER tries to load the cascading spectra
%   for ionizing collisions from the E_cascadings directory matching
%   the requested energy-array in the AURORA_root_directory. If no
%   file with matching energy-array are found new calculations are
%   done with O_e_2nd_dist, O2_e_2nd_dist and N2_e_2nd_dist, and
%   are saved in the E_cascading directory. This way the
%   unreasonable slow  calculations of these cascading-matrices
%   will only be done once for each Energy-grid configuration and
%   subsequent runs will utilize the correct pre-calculated matrices.
% 
% Calling:
%   [S2ndO,S2ndO2,S2ndN2] = e_cascading_result_finder(theta_lims2do,AURORA_root_directory)
% Input:
%  E_2do - energy-grid array (eV) [1 x (nE)]
%          double array
%  AURORA_root_directory - full path to the AURORA-root directory,
%                  automatically set with the add_AURORA script
% Output:
% 
% SEE ALSO: O_e_2nd_dist, O2_e_2nd_dist, N2_e_2nd_dist

%  Copyright © Bjorn Gustavsson 20191122, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later




% All pitch-angle-cascading files with possible matching theta_lims:
e_2nd_s_files = dir(fullfile(AURORA_root_directory,'E_cascadings'));
foundem = 0;
char(e_2nd_s_files(:).name)

for i1 = 1:numel(e_2nd_s_files),
  if ~ e_2nd_s_files(i1).isdir
    try
      load(fullfile(AURORA_root_directory,...
                    'E_cascadings',...
                    e_2nd_s_files(i1).name),...
           'E')
      if isequal(E,E_2do)
        % then we have found a match, so load
        fprintf('Loading cascading-matrices from file: %s\n',e_2nd_s_files(i1).name)
        load(fullfile(AURORA_root_directory,...
                      'E_cascadings',...
                      e_2nd_s_files(i1).name),...
             'S2ndO',...
             'S2ndO2',...
             'S2ndN2')
        foundem = 1;
        % And break the loop already
        break
      end
    catch
      % do nothing, bad practise perhaps, but I cannot be bovvered
      % to write test-exceptionings for this task of attempting to
      % load mat-files in a directory and check if the
      % possible/probable variable theta_lims matches theta_lims2do.
      % It is simple enough to just use try and see if we get a
      % match, then if there is none we still know that some work
      % needs to be done below...
    end
  end
end

if ~foundem
  % when there is no match we have to calculate the
  % cascading-matrices, there is no way around that...
  fprintf('Could not find file with matching energy grid\n')
  fprintf('Starting to calculate the requested cascading-matrices\n')
  fprintf('this will take some time.\n')
  S2ndO = O_e_2nd_dist(E,E(end),O_levels(end,1),'c');
  S2ndO2 = O2_e_2nd_dist(E,E(end),O2_levels(end,1),'c');
  S2ndN2 = N2_e_2nd_dist(E,E(end),N2_levels(end,1),'c');
  
  % ...and save the results for future use
  save_filename = sprintf('e_2nd_s_%02d_streams_%s.mat',...
                          numel(theta_lims2do)-1,...
                          datestr(now,'yyyymmdd-HHMMSS'));
  save(fullfile(AURORA_root_directory,...
                'E_cascadings',...
                save_filename),...
       'S2ndO',...
       'S2ndO2',...
       'S2ndN2')
end
