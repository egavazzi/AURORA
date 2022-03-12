function Xs = photon_ionN2(wl)
% photon_ionO - photon ionization cross section on molecular Nitrogen.
% 
% Calling:
%  Xs = photon_ionN2(wl)
% Input:
%  wl - wavelengths [1 x n] (m)
% Output:
%  Xs - ionization cross section [1 x n] (m^2)


%            Wavelength   log10(xS)
%              (A)         (cm^2)
phot_ion_N2 = [48.172      -18.758
               75.177      -18.288
               122.35      -17.738
               176.16      -17.387
               235.72      -17.084
               250.13      -17.014
               273.16      -16.968
               293.3       -16.947
               315.38      -16.895
               353.76      -16.807
               391.19      -16.727
               421.9       -16.672
               441.08      -16.654
               472.73      -16.654
               484.23      -16.657
               503.42      -16.642
               516.85      -16.629
               537.95      -16.605
               554.25      -16.623
               579.17      -16.651
               630         -16.642
               652.06      -16.639
               677.96      -16.614
               705.78      -16.59
               715.37      -16.575
               727.81      -16.66
               740.25      -16.758 ];

Xs = 10.^interp1(phot_ion_N2(:,1)*1e-10,phot_ion_N2(:,2)-4,wl,'pchip',-inf);
