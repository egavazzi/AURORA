function [t,y] = z_of_t_mirror_accup2(v0,z0,z_source,mu_avg,B,zB)

m_e = 9.1093835611e-31;       
tstart = 0;
tfinal = 0.45;% *0.22/mu_avg;
y0 = [z0; -v0*mu_avg];
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
% [t,y,te,ye,ie] = ode23(@(t,z) f2(t,z,Mu0/m_e,dBdz,zB),linspace(tstart,tfinal,250),y0,options);
[t,y] = ode23(@(t,z) f2(t,z,Mu0/m_e,dBdz,zB),linspace(tstart,tfinal,400),y0);
% Accumulate output.  This could be passed out as output arguments.

v_perp = v_perp0*sqrt(B_source/B0);
% v_par = ye(end);
% theta = atan2(v_perp,v_par);
% t_delay = te;
% mu_source = cos(theta);
% --------------------------------------------------------------------------

function dzdt = f(t,z,Mu0,B0,z0)

Re = 6.378e6;
dzdt = [z(2); -Mu0*(-(3*B0*(Re + z0)^3)/(Re + z(1))^4)];

function dzdt = f2(t,z,Mu0,dBdz,zB)

Re = 6.378e6;
dzdt = [z(2); -Mu0*interp1(zB,dBdz,z(1),'linear','extrap')];

% --------------------------------------------------------------------------

function [value,isterminal,direction] = events(t,z,z_source)
% Locate the time when z passes through z_source in a increasing direction
% and stop integration.  
value = z(1)-z_source;     % detect height = 0
isterminal = 1;   % stop the integration
direction = 1;    % positive direction
