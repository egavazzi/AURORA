function Xs = photon_ionO2(wl)
% photon_ionO - photon ionization cross section on molecular Oxygen.
% 
% Calling:
%  Xs = photon_ionO2(wl)
% Input:
%  wl - wavelengths [1 x n] (m)
% Output:
%  Xs - ionization cross section [1 x n] (m^2)


%            Wavelength   log10(xS)
%              (nm)         (cm^2)
phot_ion_O2 = [3.5698      -18.526
               7.4646      -17.952
               10.869      -17.561
               16.162      -17.202
               24.027      -16.914
               27.141      -16.832
               29.503      -16.786
               34.33      -16.755
               35.291      -16.753
               36.097      -16.746
               45.692       -16.66
               54.307      -16.596
               55.672      -16.589
               57.638      -16.624
               58.916      -16.712
               59.226      -16.722
               59.238      -16.709
               59.742      -16.598
               61.002      -16.533
               61.73      -16.523
               64.383      -16.594
               68.431      -16.626
               70.959      -16.658
               72.137      -16.596
               73.415      -16.597
               74.799      -16.827
               76.613      -17.024
               76.846      -17.032
               77.013      -17.024
               78.955      -16.914
               80.069      -16.834
               80.462      -16.842
               83.857      -17.058
               89.213      -17.401];

Xs = 10.^interp1(phot_ion_O2(:,1)*1e-9,phot_ion_O2(:,2)-4,wl,'pchip',-inf);
