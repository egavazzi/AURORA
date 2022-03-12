function Ie_omni = Ie_beams2omnidirectional(h_atm,Ie_zE)
% EXC_Z_OF_IE_ZTE - steady-state volume excitation rate
% from altitude-energy varying electron fluxes, height-varying
% neutral densities and excitation cross-section.
%   
% Calling:
%   Q = exc_z_of_Ie_zE(h_atm,E,Ie_zE,n,Xs)
% Input:
%  h_atm - height (m), double array [n_z x 1]
%  E     - energy (eV), double array [1 x n_E]
%  Ie_zE - electron particle flux (/m^2/s/Delta_E), double array
%          [n_z*n_beams x n_E]
%  n     - density of exciteable atmospheric constituent (/m^3),
%          dobule array [n_z x 1]
%  Xs    - excitation cross-section (/m^2), double array [1 x n_E]
% Output:
%  Q     - volume excitation rate (/m^3/s), time-dependent
%          variation of volume excitation-rate, double array
%          [n_z x 1]
% Example:
%  TBD
% This function calculates the excitation rates from
% electron-fluxes produced by the steady-state multistream
% electron transport-code. This function automatically determines
% the number of streams used from the size-ratio
% size(Ie_tzE,1)/size(h,1) - since the multi-stream code stacks the
% fluxes in all streams in the first dimension.

% Copyright © B. Gustavsson 20190125
%  This is free software, licensed under GNU GPL version 2 or later

szIezE = size(Ie_zE);

n_beams = size(Ie_zE,1)/numel(h_atm);
Ie_omni = 0;

if numel(szIezE) < 3
  for iB = 1:n_beams
    Ie_omni = Ie_omni + Ie_zE((1:numel(h_atm))+(iB-1)*numel(h_atm),:);
  end
else
  error('too many dimensions in Ie_zE')
end
