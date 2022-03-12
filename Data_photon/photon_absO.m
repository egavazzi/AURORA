function Xs = photon_absO(wl)
% photon_absO - photon absorption cross section on atomic Oxygen.
%  For wavelengths between 50 and 1050 A.
% Calling:
%  Xs = photon_absO2(wl)
% Input:
%  wl - wavelengths [1 x n] (m)
% Output:
%  Xs - absorption cross section [1 x n] (m^2)



%              Wavelength   ionXs/absXs
%                (m)              ()

ratio_ion2abs = [                ];


Xs = photon_ionO(wl);
