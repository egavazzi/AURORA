function ThatWentOK = Ie_zE2Qz(nN2,nO2,nO,h_atm)
% IE_ZE2Q_Z - Calculate the volume emission/excitation rates
%   Ie_zE2Q_zt loads all electron-fluxes in a directory and
%   calculates volume emission and excitation-rates for the auroral
%   emissions at 4278 Å, 6730 Å, 7774 Å, 8446 Å, and the O1D, O1S
%   and the N2A3 excitation rates. The function automatically
%   handles the calculations of excitation cross-sections etc.
% 
% Calling:
%  ThatWentOK = Ie_zE2Q_zt(nN2,nO2,nO,h_atm)
% Input:
%  nN2 - molecular nitrogen number-density (/m^3) profiles [nZ x 1]
%  nO2 - molecular oxygen number-density (/m^3) profiles [nZ x 1]
%  nO  - atomic oxygen number-density (/m^3) profiles [nZ x 1]
%  h_atm - altitude profils [m], [nZ x 1]
% Output
%  ThatWentOK - bolean output, returns 1 if processing worked out
%   

%   Copyright © 2018-2019 Bjorn Gustavsson, <bjorn.gustavsson@uit.no>
%   This is free software, licensed under GNU GPL version 2 or later

%try
  ops.display = 'quiet';
  [h_atm,E,mu_lims,Ie_ZE] = Ie_zE_loader({'.'},ops);
  szIzE = size(Ie_ZE);
  % load neutral_atm.mat
  
  emXS4278 = exc_4278(E);
  emXS6730 = exc_6730_N2(E);
  emXS8446_O = exc_8446_O(E);
  emXS8446_O2 = exc_8446_O2(E);
  emX7774_O = exc_7774_O(E);
  emX7774_O2 = exc_7774_O2(E);
  excXS_O1D = exc_O1D(E);
  excXS_O1S = exc_O1S(E);
  dE = diff(E);
  dE = dE([1:end,end]);
  XsO  = get_all_xs('O',E+dE/2);
  XsO2 = get_all_xs('O2',E+dE/2);
  XsN2 = get_all_xs('N2',E+dE/2);
  load N2_levels.dat
  load O2_levels.dat
  load O_levels.dat
  XsOi = O_levels(:,2)'*XsO;
  XsO2i = O2_levels(:,2)'*XsO2;
  XsN2i = N2_levels(:,2)'*XsN2;
  XsN2A3 = XsN2(13,:);
  
  Q4278 = exc_z_of_Ie_zE(h_atm,E(1:szIzE(2)),Ie_ZE,nN2,emXS4278(1:szIzE(2)));
  
  Q6730 = exc_z_of_Ie_zE(h_atm,E(1:szIzE(2)),Ie_ZE,nN2,emXS6730(1:szIzE(2)));
  
  Q8446 = exc_z_of_Ie_zE(h_atm,E(1:szIzE(2)),Ie_ZE,nO,emXS8446_O(1:szIzE(2))) + ...
          exc_z_of_Ie_zE(h_atm,E(1:szIzE(2)),Ie_ZE,nO2,emXS8446_O2(1:szIzE(2)));
  Q7774 = exc_z_of_Ie_zE(h_atm,E(1:szIzE(2)),Ie_ZE,nO,emX7774_O(1:szIzE(2))) + ...
          exc_z_of_Ie_zE(h_atm,E(1:szIzE(2)),Ie_ZE,nO2,emX7774_O2(1:szIzE(2)));
  
  Q8446_O = exc_z_of_Ie_zE(h_atm,E(1:szIzE(2)),Ie_ZE,nO,emXS8446_O(1:szIzE(2)));
  Q7774_O = exc_z_of_Ie_zE(h_atm,E(1:szIzE(2)),Ie_ZE,nO,emX7774_O(1:szIzE(2)));

  Q8446_O2 = exc_z_of_Ie_zE(h_atm,E(1:szIzE(2)),Ie_ZE,nO2,emXS8446_O2(1:szIzE(2)));
  Q7774_O2 = exc_z_of_Ie_zE(h_atm,E(1:szIzE(2)),Ie_ZE,nO2,emX7774_O2(1:szIzE(2)));
  
  QO1D = exc_z_of_Ie_zE(h_atm,E(1:szIzE(2)),Ie_ZE,nO,excXS_O1D(1:szIzE(2)));
  QO1S = exc_z_of_Ie_zE(h_atm,E(1:szIzE(2)),Ie_ZE,nO,excXS_O1S(1:szIzE(2)));

  QN2i = exc_z_of_Ie_zE(h_atm,E(1:szIzE(2)),Ie_ZE,nN2,XsN2i(1:szIzE(2)));
  QO2i = exc_z_of_Ie_zE(h_atm,E(1:szIzE(2)),Ie_ZE,nO2,XsO2i(1:szIzE(2)));
  QOi = exc_z_of_Ie_zE(h_atm,E(1:szIzE(2)),Ie_ZE,nO,XsOi(1:szIzE(2)));
  
  QN2A3 = exc_z_of_Ie_zE(h_atm,E(1:szIzE(2)),Ie_ZE,nN2,XsN2A3(1:szIzE(2)));
  
  save Qz_all_L.mat h_atm Q4278 Q6730 Q8446 Q7774 QO1D QO1S QN2i QN2A3 QO2i QOi Q8446_O Q7774_O Q8446_O2 Q7774_O2
  
  ThatWentOK = 1;

% catch
%  ThatWentOK = 0;
%end
