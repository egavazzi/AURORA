function Xs = photon_absO2(wl)
% photon_absO2 - photon absorption cross section on molecular Oxygen.
%  For wavelengths between 50 and 1050 A.
% Calling:
%  Xs = photon_absO2(wl)
% Input:
%  wl - wavelengths [1 x n] (m)
% Output:
%  Xs - absorption cross section [1 x n] (m^2)



%              Wavelength   ionXs/absXs
%                (m)              ()

ratio_ion2abs = [ 50e-10  0.92615
                 100e-10  0.99216
                 150e-10  0.96472
                 200e-10  1.0%085
                 250e-10  1.0%389
                 300e-10  0.98811
                 350e-10  1.0%208
                 400e-10  1.0%114
                 450e-10  0.98041
                 500e-10  1.00%37
                 550e-10  0.89249
                 600e-10  1.0%693
                 650e-10  1.0%931
                 700e-10  0.78622
                 750e-10  0.51589
                 800e-10  0.66636
                 850e-10  0.54346
                 900e-10  0.54346
                 950e-10  0];


X1 = photon_ionO2(wl(wl<=892e-10)).*interp1(ratio_ion2abs(:,1),1./ratio_ion2abs(:,2),wl(wl<=892e-10),'linear',1);
X2 = interp1([892,   949,   950,   999,   1000, 1050]*1e-10,...
             [12.817,12.817,21.108,21.108,1.346,1.346]*1e-22,...
             wl(wl>892e-10),'nearest',0);

Xs = [X1,X2];
