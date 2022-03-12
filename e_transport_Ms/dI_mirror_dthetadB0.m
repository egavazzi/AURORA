function [theta_p] = dI_mirror_dthetadB0(theta_lims,B)
% function [A_mirror,B_mirror] = dI_mirror_dthetadB(theta_lims,B)
% DI_MIRROR_DTHETADB - 
%   


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
for i_B = 1:(numel(B)-1),
  for i_th = numel(theta_lims):-1:1,
    if ~isempty(i_thetaperp) && i_th == i_thetaperp
      if B(i_B+1)/B(i_B) < 1
        theta_p(i_B,i_th) = asin(sqrt(B(i_B+1)/B(i_B))*sin(theta_lims(i_th)));
      else
        theta_p(i_B,i_th) = pi - asin(sqrt(B(i_B)/B(i_B+1))*sin(theta_lims(i_th)));
      end
    elseif theta_lims(i_th) < pi/2 % i_thetaperp
      %if B(i_B+1)/B(i_B) < 1
        theta_p(i_B,i_th) = asin(sqrt(B(i_B+1)/B(i_B))*sin(theta_lims(i_th)));
      %else
      %  theta_p(i_B,i_th) = asin(sqrt(B(i_B)/B(i_B+1))*sin(theta_lims(i_th)));
      %end
    else % theta_lims(i_th) > pi/2 % i_thetaperp
      %if B(i_B+1)/B(i_B) < 1
        theta_p(i_B,i_th) = pi-asin(sqrt(B(i_B)/B(i_B+1))*sin(theta_lims(i_th)));
      %else
      %  theta_p(i_B,i_th) = pi-asin(sqrt(B(i_B+1)/B(i_B))*sin(theta_lims(i_th)));
      %end
    end
  end
end

