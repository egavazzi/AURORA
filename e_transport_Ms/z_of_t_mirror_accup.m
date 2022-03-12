function [t_o_flight,mu_source] = z_of_t_mirror_accup(v0,z0,z_source,mu_avg,B,zB)
% z_of_t_mirror_accup - time-of-flight and pitch-angle change
% due to magnetic mirror-force.
% z_of_t_mirror_accup - calculates the time-of-flight and
% pitch-angle-cosine at distance Z1 due to magnetic
% mirror-force for electrons with initial speed V0 and
% pitch-angle-cosine MU0 at a distance Z0 when the magnetic
% field-strength are B at distance ZB.
% 
% Calling:
%  [t_o_flight,mu_source] = z_of_t_mirror_accup(v0,z0,z1,mu0,B,zB)
% Input:
%  v0  - speed (m/s), double scalar
%  z0  - distance (m), double scalar
%  z1  - end-distance (m) double scalar
%  mu0 - initial pitch-angle-cosine, double scalar
%  B   - magnetic field strength (T), double array, [nB x 1]
%  zB  - distance (m) along B, same size as B
% Output:
%  t_o_flight - time-of-flight (s), double scalar
%  mu_source  - pitch-angle-cosine at z1, double scalar
%  
% Example:
%  zB = (600:3000)*1e3;         % distance along B
%  Re = 6370e3;                 % "Earth" radius
%  B  = 5e-4*(Re./(Re+zB)).^2;  % B-field
%  E0 = 1e3;                    % Electron-energy (eV)
%  v  = v_of_E(E0)              % Electron-velocity
%  for i_mu = 1:90,
%    mu0 = cos(pi/180*(i_mu-1));
%    [t_o_flight(i_mu),mu_source(i_mu)] = z_of_t_mirror_accup(v0,...
%                                                          zB(1),...
%                                                          zB(end),...
%                                                          mu0,...
%                                                          B,...
%                                                          zB);
%  end
%  subplot(2,1,1)
%  plot(0:89,t_o_flight)
%  xlabel('Pitch-angle at z0')
%  ylabel('time (s)')
%  title('Time-of-flight between 3000 and 600 km')
%  subplot(2,1,2)
%  plot(0:89,180/pi*acos(mu_source))
%  xlabel('Pitch-angle at z0')
%  ylabel('Pitch-angle at z1')


m_e = 9.1093835611e-31;       
tstart = 0;
tfinal = 100;
y0 = [z0; v0*mu_avg];
B0 = interp1(zB,B,z0,'pchip','extrap');
B_source = interp1(zB,B,z_source,'pchip','extrap');
v_perp0 = v0*(1-mu_avg^2)^.5;
Mu0 = m_e*v_perp0^2/2/B0;

dBdz = gradient(B,zB);
refine = 4;
options = odeset('Events',@(t,z) events(t,z,z_source),'OutputSel',1,...
                 'Refine',refine);


tout = tstart;
yout = y0.';
teout = [];
yeout = [];
ieout = [];
%[t,y,te,ye,ie] = ode23(@(t,z) f(t,z,Mu0/m_e,B0,z0),[tstart tfinal],y0,options);
[t,y,te,ye,ie] = ode23(@(t,z) f2(t,z,Mu0/m_e,dBdz,zB),[tstart tfinal],y0,options);
% Accumulate output.  This could be passed out as output arguments.

v_perp = v_perp0*sqrt(B_source/B0);
v_par = ye(end);
theta = atan2(v_perp,v_par);
t_o_flight = te;
mu_source = cos(theta);
% --------------------------------------------------------------------------

function dzdt = f(t,z,Mu0,B0,z0)

Re = 6.378e6;
dzdt = [z(2); -Mu0*(-(3*B0*(Re + z0)^3)/(Re + z(1))^4)];

function dzdt = f2(t,z,Mu0,dBdz,zB)

Re = 6.378e6;
dzdt = [z(2); -Mu0*interp1(zB,dBdz,z(1),'pchip','extrap')];

% --------------------------------------------------------------------------

function [value,isterminal,direction] = events(t,z,z_source)
% Locate the time when z passes through z_source in a increasing direction
% and stop integration.  
value = z(1)-z_source;     % detect height = 0
isterminal = 1;   % stop the integration
direction = 1;    % positive direction
