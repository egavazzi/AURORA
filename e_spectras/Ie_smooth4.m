function I_of_t = Ie_smooth4(t,varagin)
% IE_SMOOTH4 - electron spectra smoothly varying in time
% with constant electron flux in all energy bins. This is for
% testing and verification of suitable time-domain behaviour of
% time-dependent electron transpotr functions.
%   
% 
% Calling:
%  Ie_of_t = Ie_smooth4(t[,E])
% Input:
%  t - time (s), double array [1 x nt], output fluxes are well
%      defined between 0 and 0.3 s.
%  E - energy grid, optional unused input argument
% Output:
%  Ie - differential electron flux (#/m^2/s)
% 

% Copyright © B Gustavsson 20180502, bjorn.gustavsson@uit.no
% This is free software, licensed under GNU GPL version 2 or later


t0 = [0,     2,     29,    39,    49,    69,    79,    99,   199,300]/1000;
I0 = [5e+10, 1e+11, 1e+11, 3e+11, 2e+11, 2e+11, 1e+11, 5e+11,  0,  0];


I_of_t = interp1(t0,I0,t,'pchip');