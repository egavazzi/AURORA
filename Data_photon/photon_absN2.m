function Xs = photon_absN2(wl)
% photon_ionO - photon absorption cross section on molecular Nitrogen.
% 
% Calling:
%  Xs = photon_absN2(wl)
% Input:
%  wl - wavelengths [1 x n] (m)
% Output:
%  Xs - ionization cross section [1 x n] (m^2)



%              Wavelength   ionXs/absXs
%                (m)              ()

ratio_ion2abs = [ 50e-10    0.77063
                 100e-10    0.85959
                 150e-10    0.82641
                 200e-10    0.88314
                 250e-10    1.0
                 300e-10    0.96937
                 350e-10    1.0
                 400e-10    0.98155
                 450e-10    0.94581
                 500e-10    0.97572
                 550e-10    1.0
                 600e-10    0.97458
                 650e-10    0.75874
                 700e-10    0.76855
                 750e-10    0.76855
                 800e-10    0
                 850e-10    0
                 900e-10    0
                 950e-10    0];

X1 = photon_ionN2(wl(wl<=740e-10)).*interp1(ratio_ion2abs(:,1),1./ratio_ion2abs(:,2),wl(wl<=740e-10),'linear',1);
X2 = interp1([700,   749,   750,   799,   800,   849,   850,   899,   900, 949, 950,   1000]*1e-10,...
             [24.662,24.662,33.578,33.578,16.992,16.992,20.249,20.249,9.68,9.68,50.988,50.988]*1e-22,...
             wl(wl>740e-10),'nearest',0);

Xs = [X1,X2];
