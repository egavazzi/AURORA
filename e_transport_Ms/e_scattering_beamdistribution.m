function [Pmu2mup,theta2beamW,BeamW] = e_scattering_beamdistribution(mu_lims,n_dirs)
% e_scattering_beamdistribution - angular redistribution PDF for
% scattering from all pitch-angles to discrete beams with arbitrary
% pitch-angle limits, and the relative contribution-weighting-
% matrices accounting for isotropic distribution of fluxes inside
% corresponding beams.
% 
% Calling:
%  [Pmu2mup,theta2beamW] = e_scattering_beamdistribution(mu_lims,n_dirs)
% Input:
%   mu_lims - cosine of pitch-angle limits, double array [1 x (nB-1)]
%             or [(nB-1) x 1 ], -1 and +1 are the effective min and
%             maximum cosine of pitch-angle limits
%   n_dirs  - number of directions to use for the discretized
%             calculation, optional argument, when left emtpy
%             defaults to 181
% Output:
%  Pmu2mup     - double array [n_dirs x n_dirs x n_beams] with the
%                probability distribution for an electron ending up
%                in a beam starting at pitch-angle theta0 and
%                scattering an angle theta1.
%  theta2beamW - Relative weighting matrix, double matrix [n_beams
%                x n_dirs], with the relative contribution from
%                within each beam, assuming isotropic pich-angle
%                distribution within each beam.
%  BeamW       - pitch-angle-stream solid angle array (ster) double
%                array [1 x n_mu]
% Example:
%   n_dirs = 181;
%   mu_lims = [-1 -cos(pi/4) 0 cos(pi/4) 1];
%   [Pmu2mup,theta2beamW,BeamW] = e_scattering_beamdistribution(mu_lims,n_dirs);
%   for i1 = 1:4
%     subplot(3,2,i1)
%     imagesc(0:180,0:180,Pmu2mup(:,:,i1))
%   end
%   subplot(3,1,3)
%   imagesc(0:180,1:(numel(mu_lims+1)),theta2beamW)
%   

%  Copyright � Bjorn Gustavsson 20180603, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later



if nargin < 2 || isempty(n_dirs)
  n_dirs = 181;
end

theta0 = linspace(0,pi,n_dirs); % Initial pitch-angle
theta1 = linspace(0,pi,n_dirs); % scattering angle, i.e. the angle
                                % a particle scatters of its
                                % initial traveling direction,
                                % without taking the direction
                                % around that initial direction
                                % into account
theta0 = mu_avg(0:(180/(n_dirs-1)):180)*pi/180;
theta1 = mu_avg(0:(180/(n_dirs-1)):180)*pi/180;

for i1 = numel(theta0):-1:1
  if i1 == numel(theta0) || mod(i1,10) == 0
    fprintf('starting with %d/(of %d) at: %s\n',i1,numel(theta0),datestr(now,'HH:MM:SS'))
  end
  % First we calculate the unit vector of initial propagation
  e0 = [sin(theta0(i1)) 0 cos(theta0(i1))]';
  
  for i2 = numel(theta1):-1:1
    % Then we rotate that one with a(ll) scattering angles:
    e1 = rot_around_v([0 1 0],theta1(i2))*e0;
    for iPhi = 1:(n_dirs-1)
      % Then we rotate that one around all clock-angles
      es = rot_around_v(e0',180/(n_dirs-1)*iPhi*pi/180)*e1;
      % then we save away the corresponding final pitch-angle of
      % the unit vector after scattering:
      mu(iPhi) = es(3);
    end
    % So the number of mu (cosine of pitch-angles) that are within
    % the beam pitch-angle limits ...
    for iMu = (numel(mu_lims)-1):-1:1
      B(i1,i2,iMu) = sum(mu_lims(iMu)<=mu&mu<mu_lims(iMu+1));
    end
  end
end
% ...can be converted to the fraction - which should be equal to
% the probability P of going from theta to theta'
Pmu2mup = B./repmat(sum(B,3),[1 1 size(B,3)]);
Pmu2mup(1,1,:) = Pmu2mup(2,1,:)/2+Pmu2mup(1,2,:)/2; 
Pmu2mup(end,end,:) = Pmu2mup(end-1,end,:)/2+Pmu2mup(end,end-1,:)/2; 
% B(1,from_mu,to_mup)
for iMu = (numel(mu_lims)-1):-1:1
  theta2beamW(iMu,:) = abs(sin(theta0)).*(mu_lims(iMu)<cos(theta0)&cos(theta0)<=mu_lims(iMu+1));
  BeamW(iMu) = abs(integral(@(pa) sin(pa),acos(mu_lims(iMu)),acos(mu_lims(iMu+1))));
end
