function [t_delay,mu_end] = dt_of_z_mirror(v0,z0,z_end,mu_avg,B,zB)
% dt_of_z_mirror - travel-time along magnetic field with mirror-force 
% This function calculates electron-travel-time and change in 
% pitch-angle-cosine between star-altitude and a
% end-altitude. (Typically from the top of the ionosphere to a
% source-altitude)
% 
% Calling:
%  [t_delay,mu_end] = dt_of_z_mirror(v0,z0,z_end,mu_avg,B,zB)
% Input:
%  v0     - initial speed (m/s), double scalar
%  z0     - initial altitude (m), double scalar
%  z_end  - end-altitude (m), double scalar
%  mu_avg - initial pitch-angle-cosine, double scalar
%  B      - magnetic-field strength (T), double array [nB x 1]
%  zB     - array with altitudes (m) for B, double array [nB x 1]
% Output:
%  t_delay - travel-time (s) from z0 to z_end, double scalar
%  mu_end  - pitch-angle-cosine at z_end, double scalar
%
% Example:
%  zB = (600:3000)*1e3;         % distance along B
%  Re = 6370e3;                 % "Earth" radius
%  B  = 5e-4*(Re./(Re+zB)).^2;  % B-field
%  E0 = 1e3;                    % Electron-energy (eV)
%  v  = v_of_E(E0)              % Electron-velocity (m/s)
%  for i_mu = 1:90,
%    mu0 = cos(pi/180*(i_mu-1));
%    [t_o_flight(i_mu),mu_source(i_mu)] = dt_of_z_mirror(v_of_E(E_0),z0,z_end,cos((i1-1)*pi/180),B0,lRB);
%  end
%  subplot(2,1,1)
%  plot(0:89,t_o_flight)
%  xlabel('Pitch-angle at z0')
%  ylabel('time (s)')
%  title('Time-of-flight between 3000 and 600 km')
%  subplot(2,1,2)
%  plot(0:89,180/pi*acos(mu_source))
%  xlabel('Pitch-angle at z0')
%  ylabel('Pitch-angle at z_end')

%   Copyright © 2020 Bjorn Gustavsson, <bjorn.gustavsson@uit.no>
%   This is free software, licensed under GNU GPL version 2 or later

m_e = 9.1093835611e-31;    % Electron-mass (kg)
tstart = 0;                % Start-time (s)
tfinal = 1000;             % Maximum end-time (s)

% Magnetic-field at start and end-points, interpolate to improve accuracy
B0 = interp1(zB,B,z0,'pchip','extrap');
B_end = interp1(zB,B,z_end,'pchip','extrap');

% Magnetic-field-strength-gradient
dBdz = gradient(B,zB);

% Initial conditions, altitude and field-aligned velocity-component
y0 = [z0; v0*mu_avg];

% Initial perpendicular velocity and initial magnetic-moment
v_perp0 = v0*(1-mu_avg^2)^.5;
Mu0 = m_e*v_perp0^2/2/B0;

% Options for ode23, primarily the event for passing goal-line
refine = 4;
options = odeset('Events',@(t,z) events(t,z,z_end),'OutputSel',1,...
                 'Refine',refine);

% tout = tstart;
% yout = y0.';
% teout = [];
% yeout = [];
% ieout = [];

% [t,y,te,ye,ie] = ode23(@(t,z) f2(t,z,Mu0/m_e,dBdz,zB),[tstart tfinal],y0,options);
[~,~,te,ye] = ode23(@(t,z) f2(t,z,Mu0/m_e,dBdz,zB),[tstart tfinal],y0,options);

% Accumulate output.
v_perp = v_perp0*sqrt(B_end/B0);
v_par = ye(end);
theta = atan2(v_perp,v_par);
t_delay = te;
mu_end = cos(theta);


% --------------------------------------------------------------------------
% $$$ function dzdt = f(~,z,Mu0,B0,z0)
% $$$ 
% $$$ Re = 6.378e6;
% $$$ dzdt = [z(2); -Mu0*(-(3*B0*(Re + z0)^3)/(Re + z(1))^4)];

% function dzdt = f2(t,z,Mu0,dBdz,zB)
function dzdt = f2(~,z,Mu0,dBdz,zB)

dzdt = [z(2); -Mu0*interp1(zB,dBdz,z(1),'pchip','extrap')];

% --------------------------------------------------------------------------

function [value,isterminal,direction] = events(~,z,z_end)
% Locate the time when z passes through z_end in a increasing direction
% and stop integration.  
value = z(1)-z_end;     % detect height = 0
isterminal = 1;   % stop the integration
direction = 1;    % positive direction
