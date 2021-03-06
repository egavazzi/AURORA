function sigma = exc_8446_O2(E)
% Xs = exc_8446_O2(E)
% 
% exc_8446 - emission cross section (m^2) for 8446 AA
% from molecular oxygen. E electron energy (eV)
%
% Electron-impact production of OI 3p 3P 8446 A transition
% in the dissociative ionization/excitation of O2
% Schulman et al, Phys. Rev. A, 32, 2100, 1985
% 
% Digitised from fig.7 by DKW March 2009. Value in paper is 2 +/- 15% @ 100eV.
% parent       O2
% products     0    (emission)
% threshold    15.9   (???) - seems low, dissociation energy of O2
%                             is 5.16 eV and the energy-level of
%                             O(3p3P) is 10.99 eV, so threshold has
%                             to  be 16.15 eV /BG 20180529
% units        eV
% valid from   0.0   Remove this eventually - or replace with other comments ...
% valid to     0.0
% points       39

EnX = [ 15.9877      0.00888889
        16.3801      0.124444
        16.7724      0.24
        17.1341      0.377778
        18.1027      0.475556
        19.0345      0.6
        21.5111      0.804444
        24.0429      0.968889
        26.5992      1.11556
        31.3502      1.27111
        34.5563      1.34667
        37.787       1.40444
        41.5816      1.45333
        42.6299      1.49333
        44.7387      1.56444
        48.4843      1.64889
        51.1387      1.72444
        57.6123      1.83111
        64.6682      1.91556
        70.6697      1.96444
        77.2536      1.99111
        84.9594      2.00444
        93.2475      1.99556
        103.767      1.96889
        114.863      1.92444
        130.973      1.84444
        147.09       1.76
        159.877      1.68889
        178.207      1.6
        195.96       1.52889
        215.375      1.45333
        231.448      1.4
        253.051      1.33778
        279.08       1.26667
        301.769      1.21778
        327.221      1.16444
        348.794      1.12444
        375.344      1.07556
        398.572      1.03556
        1e3          0.58735   % The last 3 rows added by B Gustavsson
        1e4          0.14193   % 201806 25 to make extrapolation to 
        1e5          0.034298];% higher energies well-behaved
% These values are based on polynomial etrapolation of log of cross-section
% vs log of energy

sigma = exp(interp1(log(EnX(:,1)),log(EnX(:,2)),log(E),'pchip'));
sigma(E<16.15) = 0;
sigma = sigma*1e-22;
