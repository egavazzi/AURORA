function [Mi2ip1,Mi2im1,theta_p] = dI_mirror_dthetadB(theta_lims,B)
% DI_MIRROR_DTHETADB - Magnetic mirroring multi-stream transfer
% DI_MIRROR_DTHETADB calculates the pitch-angle-stream-to-
% pitchangle-stream transfer caused by the magnetic mirror-force. 
% The transfer between pitch-angle-streams from one altitue to
% another, with different magnetic field-strengths are calculated
% using B_2/B_1 = sin(theta_2)^2/sin(theta_1)^2 and assuming that
% the pitch-angle distribution is isotropic inside every
% pitch-angle-stream at every step
%
% Calling:
%   [Mi2ip1,Mi2im1,theta_p] = dI_mirror_dthetadB(theta_lims,B)
% Input:
%   theta_lims - pitch-angle-limits for the streams used (radians),
%                double array [1 x (n_streams+1)], typically in the
%                range [0,pi]
%   B          - magnetic-field (T), either double array [n_z x 1] 
%                or [n_z x 3], if the magnetic field is given
%                instead of the magnetic field-strength the
%                function will adapt.
% Output:
%   Mi2ip1  - cell-array [1 x n_mu] with transfer-factor-vectors
%             [n_z x 1] from pitch-angle-stream(s) i_mu to i_mu+1
%   Mi2im1  - cell-array [1 x n_mu] with transfer-factor-vectors
%             [n_z x 1] from pitch-angle-stream(s) i_mu to i_mu-1
%   theta_p - pitch-angle-limit-chage-array [n_z x (n_mu+1)]
%             containing the pitch-angle-limit at altitude i_z
%             mirrored to altitude i_z+1
% 
% See also: Ie_Mstream_tz_Bm_2_aurora_faster, Ie_M_stream_2_photo_e
%           de_M_stream_CNztusDB, de_M_stream2HS_us 

%  Copyright � Bjorn Gustavsson 20200220, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later

if size(B,2) == 3
  % Then we have sent in the magnetic-field 
  % but we need magnetic field strength
  B = ( B(:,1).^2 + B(:,2).^2 + B(:,3).^2 ).^(1/2);
end

[dtheta_p,i_thetaperp] = min(abs(theta_lims-pi/2));
if dtheta_p > 1e-5
  % then we for some reason have no theta_lims perp to B
  % TODO: make up a parallel/anti-parallel splitting scheme here
  i_thetaperp = [];
end

dB = diff(B);
% This will start to fail if the close-to-perpendicular
% pitch-angle-limits start to become narrow compared to the change
% of pitch-angle due to mirroring between adjacent magnetic-field
% strengths. 
% TODO: find a fix for this
% TODO: at least make a warning for when this happens.
% TODO: vectorize this, but no hurry just to be cute.
for i_B = 1:(numel(B)-1)
  for i_th = numel(theta_lims):-1:1
    % This is the change of pitch-angles due to magnetic mirroring
    % applied for the pitch-angle-limits of the streams, only.
    if ~isempty(i_thetaperp) && i_th == i_thetaperp
      if B(i_B+1)/B(i_B) < 1
        theta_p(i_B,i_th) = asin(sqrt(B(i_B+1)/B(i_B))*sin(theta_lims(i_th)));
      else
        theta_p(i_B,i_th) = pi - asin(sqrt(B(i_B)/B(i_B+1))*sin(theta_lims(i_th)));
      end
    elseif theta_lims(i_th) < pi/2
      theta_p(i_B,i_th) = asin(sqrt(B(i_B+1)/B(i_B))*sin(theta_lims(i_th)));
    else % theta_lims(i_th) > pi/2 % i_thetaperp
      theta_p(i_B,i_th) = pi-asin(sqrt(B(i_B)/B(i_B+1))*sin(theta_lims(i_th)));
    end
  end
end

for i_stream = 1:(numel(theta_lims)-1)
  
  % Here we divide the proceedings such that for B decreasing the
  % mirroring will lead to shift to next higher stream and to the
  % previous lower stream. For the case where 180 degrees are
  % field-aligned down on the northern hemisphere in the first
  % stream and consecutive streams are less and less downward this
  % means that Mi2ip1 is the loss coefficient for flux in stream i
  % to stream i+1. Likewise Mi2im1 is the loss coefficient for flux
  % in stream i to stream i-1.
  Mi2ip1{i_stream} = max(0, (cos(theta_lims(i_stream+1)) - cos(theta_p(:,i_stream+1)))./...
                            (cos(theta_p(:,i_stream)) - cos(theta_p(:,i_stream+1))) );
  
  Mi2im1{i_stream} = max(0, -(cos(theta_lims(i_stream)) - cos(theta_p(:,i_stream)))./...
                         (cos(theta_p(:,i_stream)) - cos(theta_p(:,i_stream+1))));
  % (int(sin(th),th,theta_ip1,thetaP_ip1))/int(sin(th),th,thetaP_i,thetaP_ip1) 
  % 
  %  (cos(theta_ip1) - cos(thetaP_ip1))/(cos(thetaP_i) - cos(thetaP_ip1))
  
end
if ~isempty(i_thetaperp)
  Mi2ip1{i_thetaperp-1} = Mi2ip1{i_thetaperp-1}/2;
  Mi2im1{i_thetaperp} = Mi2im1{i_thetaperp}/2;
end
% keyboard