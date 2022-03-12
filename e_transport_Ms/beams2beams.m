function B2B = beams2beams(phase_fcn,Pmu2mup,thetaW)
% BEAMS2BEAMS - multi-beam scattering redistribution coefficients 
%   
% Calling:
%   B2B = beams2beams(phase_fcn,Pmu2mup,thetaW)
% Input:
%  phase_fcn - phase-function, [n_dirs x 1], double array with
%              probabillity for scattering with an angle theta, for
%              scattering normally probability is given for
%              scattering to a solid angle dOhmega at an angle
%              theta, here that probabillity should be weighted
%              with 2*pi*sin(theta)*dtheta
%  Pmu2mup   - array with probabillities for scattering to beams
%              from arbitrary pitch-angles to the different beams
%              given a scattering of theta, double array
%              [n_dirs x n_dirs x nBeams], see:
%              e_scattering_beamdistribution 
%  thetaW    - Relative weighting matrix, double matrix [n_beams
%              x n_dirs], with the relative contribution from
%              within each beam, assuming isotropic pich-angle
%              distribution within each beam.
% Output:
%  B2B       - beam-redistribution coefficients.
% Example:
%   n_dirs = 181;
%   mu_lims = [-cos(pi/4) 0 cos(pi/4)];
%   theta = linspace(0,pi,181);
%   E = 2e3;
%   [Pmu2mup,theta2beamW] = e_scattering_beamdistribution(mu_lims,n_dirs);
%   for i1 = 1:4
%     subplot(3,2,i1)
%     imagesc(0:180,0:180,Pmu2mup(:,:,i1)),colorbar_labeled('')
%   end
%   subplot(3,2,1)
%   xlabel('scattering angle')
%   ylabel('starting angle')
%   title('to beam #1')
%   subplot(3,2,5)
%   imagesc(0:180,1:(numel(mu_lims+1)),theta2beamW),colorbar_labeled('')
%   [phfcnN2_e,phfcnN2_i] = phase_fcn_N2(E,theta');
%   phfcnN2_e = phfcnN2_e.*sin(theta')./sum(phfcnN2_e.*sin(theta'));
%   B2B = beams2beams(phase_fcn,Pmu2mup,thetaW)
%   subplot(3,2,6)
%   imagesc(B2B),colorbar_labeled('')

%  Copyright © Bjorn Gustavsson 20180603, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later


for i3 = size(Pmu2mup,3):-1:1,
  B2B(i3,:) = thetaW*(Pmu2mup(:,:,i3)*phase_fcn(:));
end
B2B = B2B./repmat(sum(B2B,1),size(B2B,1),1);
