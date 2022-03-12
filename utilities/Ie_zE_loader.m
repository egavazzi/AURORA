function [h_atm,E,mu_lims,IeZE] = Ie_zE_loader(data_path,ops)
% Ie_zE_loader reads flickering-aurora data from the
% time-dependent electron transport runs for flickering
% aurora. This run was done in sequences of short time-intervalls,
% typically 0.1 s, with initial fluxes from the end of the
% preceding time-period. To ease further analysis we need to strip
% off the fluxes at the first time-step from data from the second
% file onwards. This function allows loading of data from two
% separate runs, stored in separate directories for addition or
% subtraction of the fluxes - this provided that the altitude,
% pithc-angle, time and energy-grids all match up
% 
% Calling:
%  [t,h_atm,E,mu_lims,IeZE] = Ie_zE_loader(data_path)
% Input:
%  data_path - path, absolute or relative path to directory with
%              output from the time-dependent electron transport
%              runs for flickering aurora, cell array [1 x 1] or
%              [1 x 2]
% Output:
%  h_atm    - altitude (m), double array [n_z x 1]
%  E        - energy grid (eV), double array [1 x n_E]
%  mu_lims  - cosine-of-pitch-angle-limits for the number of
%             streams, double array [1 x (n_streams+1)]
%  IeZ E    - electron number flux (#/m^2/s), double array 
%             [(n_z*n_streams) x n_t x n_e] with number of
%             electrons per energy-bin, sorted along the first
%             dimension as stream-above-stream from most paralell
%             to B downward to most parallel to B upward.

%   Copyright © 20190601 Bjorn Gustavsson, <bjorn.gustavsson@uit.no>
%   This is free software, licensed under GNU GPL version 2 or later



  
%  E               1x94                     752  double
%  Ie_zE       4833x31x94            112666896  double
%  h_atm         537x1                     4296  double
%  mu_lims         1x10                      80  double
if ~strcmp(ops.display,'quiet')
  res_files = dir(fullfile(data_path{1},'Ie*.mat'))
  disp(['loading: ',fullfile(data_path{1},res_files(1).name)])
else
  res_files = dir(fullfile(data_path{1},'Ie*.mat'));
end
load(fullfile(data_path{1},res_files(1).name))
IeZE = Ie_zE;
for i1 = 2:length(res_files),
  if ~strcmp(ops.display,'quiet')
    disp(['loading: ',fullfile(data_path{1},res_files(i1).name)])
  end
  load(fullfile(data_path{1},res_files(i1).name),...
       'E','Ie_zE','h_atm','mu_lims','t_run')
  IeZE(:,1:size(Ie_zE,2),i1) = Ie_zE;
end
