function [t_delay,mu_source] = z_of_t_mirror(v0,z0,z_source,mu_avg,B,zB)
%BALLODE  Run a demo of a bouncing ball.  
%   This is an example of repeated event location, where the initial
%   conditions are changed after each terminal event.  This demo computes ten
%   bounces with calls to ODE23.  The speed of the ball is attenuated by 0.9
%   after each bounce. The trajectory is plotted using the output function
%   ODEPLOT. 
%
%   See also ODE23, ODE45, ODESET, ODEPLOT, FUNCTION_HANDLE.

%   Mark W. Reichelt and Lawrence F. Shampine, 1/3/95
%   Copyright 1984-2004 The MathWorks, Inc.
%   $Revision: 1.17.4.2 $  $Date: 2005/06/21 19:24:10 $

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
t_delay = te;
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
