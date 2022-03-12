function v = v_of_E(E)
% v_of_E - convert electron energy (eV) to velocity (m/s)
%
% Calling:
%  v = v_of_E(v)
% Input:
%  E - electron energy (eV)  double [ N x M ]
% Output:
%  v - electron velocity (m/s) double [ N x M ]


%  Copyright © Bjorn Gustavsson 20030403, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later

m_e     = 9.10939e-31;
q_e     = 1.6021773e-19;

v = (2*q_e*abs(E)/m_e).^(1/2).*sign(E);
