function D_e = time_of_flight_DOMS_diffusion(E,dE,theta_lims,plotit)
% TIME_OF_FLIGHT_DOMS_DIFFUSION - Energy-pitch-angle time-of-flight widening
%   of electrons. For discrete ordinate electron-multi-beam codes
%   the electron flux (Ie(z,E,theta,t) is discretized into energy
%   bins and pitch-angle streams. With the assumption that the
%   electron flux is isotropically distributed inside an
%   energy-pitch-angle bin (flux with energies between E(i) and
%   E(i+1) with pitch-angles to the magnetic field (B) between
%   theta(j) and theta(j+1)) the electrons will have a distribution
%   of field-aligned velocity components v(E)*cos(theta). That
%   spread in velocities will lead to a time-of-flight
%   dispersion. This function calculates the width of that t-o-f
%   dispersion and converts that to a diffusion coefficient D.
%   
% Calling:
%   D_e = time_of_flight_DOMS_diffusion(E,dE,theta_lims,plotit)
% Input:
%  E          - Energy (eV), double array [nE x 1], 
%  dE         - width in energy (eV) of energy bins, double array
%               [nE x 1] 
%  theta_lims - pitch-angle limits (radians) of streams, double
%               array [nStreams + 1 x 1]
%  plotit     - Flag for plotting the proceedings, Bolean scalar,
%               optional input, default is 0 i.e. not to plot 
% Output:
%  D_e        - Diffusion-coefficients, double matrix [nE x nStreams]

%  Copyright © Bjorn Gustavsson 20190515, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later


if nargin < 4 || isempty(plotit)
  plotit = 0;
  nE = 3;
  nTheta = 3;
  n_ti = 3;
  n_thi = 4;
else
  nE = 343;
  nTheta = 131;
  n_ti = 701;
  n_thi = 401;
end

for iE = numel(E):-1:1
  v = linspace(v_of_E(E(iE)),v_of_E(E(iE)+dE(iE)),nE);
  for iTheta = 1:(numel(theta_lims)-1)
    theta_a = theta_lims(iTheta);
    theta_b = theta_lims(iTheta+1);
    if theta_lims(iTheta) == pi/2
      % theta_a = theta_lims(iTheta)*0.9 + 0.1*theta_lims(iTheta+1);
      theta_a = theta_lims(iTheta)*0.8 + 0.2*theta_lims(iTheta+1);
    end
    if theta_lims(iTheta+1) == pi/2
      % theta_b = theta_lims(iTheta)*0.1 + 0.9*theta_lims(iTheta+1);
      theta_b = theta_lims(iTheta)*0.2 + 0.8*theta_lims(iTheta+1);
    end
    theta = linspace(theta_a,theta_b,nTheta);
    theta4i = linspace(min(theta),max(theta),n_thi);
    [V,Theta] = meshgrid(v,theta);
    v_par = V.*cos(Theta);
    t_arrival = 500e3./v_par;
    at_a = (max(t_arrival(:)) + min(t_arrival(:)))/2;
    dt_a = max(t_arrival(:)) - min(t_arrival(:));
    D = (dt_a/4)^2/at_a;
    D_e(iE,iTheta) = abs(D);
    if plotit
      W = sin(Theta);
      W(:,[1 end]) = 0;
      t4i = linspace(min(t_arrival(:)),max(t_arrival(:)),701);
      [T4i,Theta4i] = meshgrid(t4i,theta4i);
      Warrival = griddata(t_arrival(:),Theta(:),W(:),T4i,Theta4i);
      clf
      subplot(3,1,1)
      pcolor(t_arrival,Theta,W),shading flat   
      subplot(3,1,2)
      pcolor(T4i,Theta4i,Warrival),shading flat
      subplot(3,1,3)
      plot(T4i(1,:),nansum(Warrival))
      hold on,
      plot(T4i(1,:),max(nansum(Warrival))*exp(-(T4i(1,:)-at_a).^2/(dt_a/4)^2),'r'),
      drawnow
    end
  end
end
