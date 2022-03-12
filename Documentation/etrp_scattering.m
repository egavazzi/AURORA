%% Script illustrating (elastic) scattering
%  For electron transport collisions lead to change in propagation
%  velocity of electrons, for ellastic scattering direction of
%  propagation is modified (for inelastic scattering also the speed
%  since kinetic energy is lost to excitation/ionization). This
%  script illustrates how the angular change of
%  propagation-direction leads to change of pitch-angle

%%  Physical constants
m_e     = 9.1093835611e-31;       % electron rest mass [kg]
q_e     = 1.602176620898e-19;     % elementary charge [C]

%% Magnetic field and electron gyro-frequency
%  of typical ionospheric strength
B = 5e-5; % T
w_e = w_e_gyro(B); % electron angular gyro-frequency (rad/s)
f = w_e/2/pi;      % (Hertz)
T = 1/f;           % gyro-period (s)

%% Initial electron speed and gyro-radius:
% For illustrating purposes we can chose an electron with 10 eV
% energy
v1 = v_of_E(10); % m/s
r_g = v1/w_e;   % m

%% Initial pitch-angle and velocity-vector
% Any initial pitch-angle will do, a reasonably large pitch-angle
% makes the illustration clear
theta = 80/180*pi;
v0 = v1*[sin(theta),0,-cos(theta)];

%% Scattering angle
psi = -50/180*pi;

%% Pitch-angle discretization into finite number of streams:
if ~exist('mu_lims','var')
  mu_lims = -1:(1/4):1;
end
%% Particle trajectory:
%  calculated with the Boris-mover scheme
rv = ode_boris_mover([0:250]/100*T,[0,r_g,1,v0],q_e,m_e,[0 0 0]',B*[0 0 1]',w_e);

%% Velocity-vector to plot, suitable scaling to magnitude similar to
% gyro-radius
v2p = rv(4:6,end)/norm(v0)/3;
%% Velocity-change due to scattering
vp = (rot_around_v([0 0 1],psi)*v2p);
for i1 = 1:360,
  v4psi(i1,:) =  (rot_around_v(v2p/norm(v2p),i1/180*pi)*vp);
end
mu = v4psi(:,3)/norm(vp);
subplot(1,3,1)
plot3(rv(1,:),rv(2,:),rv(3,:),'k','linewidth',2)
hold on
plot3(rv(1,end),rv(2,end),rv(3,end),'g.','markersize',22)
arrow(rv(1:3,end)',rv(1:3,end)'+v2p','width',1)
arrow(rv(1:3,end)',rv(1:3,end)'+vp',...
      'edgecolor','m','facecolor','m','width',1)

R_arc = r_arc(v2p',vp',rv(1:3,end),0.2);
plot3(R_arc(:,1),R_arc(:,2),R_arc(:,3),'k')
axis equal
xt = -0.3;
yt = -0.1;
zt = 0.4; 
text(xt,yt,zt,'\psi','fontsize',14)
axis([-0.34038      0.20     -0.46021 0.2 0.13198 1])

subplot(1,3,2)
plot3(rv(1,:),rv(2,:),rv(3,:),'k','linewidth',2)
hold on
plot3(rv(1,end),rv(2,end),rv(3,end),'g.','markersize',22)
arrow(rv(1:3,end)',rv(1:3,end)'+v2p','width',1)
arrow(rv(1:3,end)',rv(1:3,end)'+vp',...
      'edgecolor','m','facecolor','m','width',1)
hcl = cline(rv(1,end)+v4psi(:,1),...
            rv(2,end)+v4psi(:,2),...
            rv(3,end)+v4psi(:,3),...
            mu,jet,[-1,1]);
axis equal
colorbar_labeled('\mu')
axis([-0.34038      0.20     -0.46021 0.2 0.13198 1])

subplot(1,3,3)
hcl2 = cline(v4psi(:,1)./sum(v4psi.^2,2).^.5,...
             v4psi(:,2)./sum(v4psi.^2,2).^.5,...
             v4psi(:,3)./sum(v4psi.^2,2).^.5,...
             mu,jet,[-1,1]);
view(90,0)
hold on

plot3([0,v2p(1)'/norm(v2p)],...
      [0,v2p(2)'/norm(v2p)],...
      [0,v2p(3)'/norm(v2p)],...
      'k+')

view(90,0)
% axis equal
plot3([-1 1],[-1 1],mu_lims'*[1 1],'k')
mu_centres = mu_lims(1:end-1)/2 + mu_lims(2:end)/2;
for i1 = numel(mu_centres):-1:1
  mu_names{i1} = sprintf('u_%d',i1);
end
axis equal
% $$$ set(gca,...
% $$$     'xtick',[],...
% $$$     'ytick',[],...
% $$$     'ztick',mu_centres,...
% $$$     'zticklabel',{'u_1','u_2','u_3','u_4','u_5','u_6','u_7','u_8'})
set(gca,...
    'xtick',[],...
    'ytick',[],...
    'ztick',mu_centres,...
    'zticklabel',mu_names)
axis([-1 1 -1 1 -1 1])
drawnow

qwe1 = mfilename('fullpath');
[qweA,qweB,qweC] = fileparts(qwe1);

print('-depsc2','-painters',fullfile(qweA,'e-scattering-geometry-01.eps'));
subplot(1,3,1)
grid on
subplot(1,3,2)
grid on
print('-depsc2','-painters',fullfile(qweA,'e-scattering-geometry-02.eps'));

vaz = -65.5;
vel =  16;
subplot(1,3,1)
view(vaz,vel);
subplot(1,3,2)
view(vaz,vel);
print('-depsc2','-painters',fullfile(qweA,'e-scattering-geometry-04.eps'));
subplot(1,3,1)
grid off
subplot(1,3,2)
grid off
print('-depsc2','-painters',fullfile(qweA,'e-scattering-geometry-03.eps'));

