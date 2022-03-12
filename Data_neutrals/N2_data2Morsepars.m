function [a,nu0] = N2_data2Morsepars(dE0to1,De,m,r_e)
% N2_DATA2MORSEPARS - 
%   

q_e     = 1.602176565e-19;  % elementary charge [C]
a0 = 3;

a = fzero(@(a) errfcn(a,De,r_e,m,dE0to1),a0 );

nu0 = a*1e10*sqrt(2*De/q_e/m)/(2*pi)*q_e;

function dE = errfcn(a,De,r_e,m,dE0to1)
% ERRFCN - 
%   

[~,E_vib] = potential_Morse(r_e,[0,1],De,r_e,a,0,m);

dE = diff(E_vib) - dE0to1;