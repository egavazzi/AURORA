function [Ie_zt] = de_M_stream_CNzt(h_atm,t_in,mu,I_top_of_t,I0,v,QC_mu,A,B_b2b,Dz,lowerBvals)
% de_M_stream_CNzt - multistream Crank-Nicholson electron transport
% solver with central difference operator for spatial
% derivatives. This integrates the time-dependent electron
% transport equation:
% 
%  1/v(E)*DI_mu/Dt = -mu*DI_mu/Dz -A*I_mu + B*I_mu'
% 
% For electron fluxes at energy E in an arbitrary number of
% pitch-angle streams with average cosine-of-pitch-angle mu. The
% altitude profiles for the initial fluxes in I0 are to be stacked
% ontop of each other, [I0_mu1;I0_mu2;...;I0_mun], and the upper
% boundary condition I_top_of_t is to be a cell-array with the
% upper boundary condition as a function of time pitch-angle by
% pitch-angle.
% 
% Calling:
%  [Ie_zt] = de_M_stream_CNzt(h_atm,t,mu,I_top_of_t,I0,v,...
%                               QC_mu,A,B,Pmu2mup,Dz,...
%                               lowerBvals)
% Input:
%   h_atm - distance along themagnetic field (tilted altitude) [nz x 1] (m)
%   t     - time, [1 x nt] (s)
%   mu    - pitch-angle [1 x n_mu] (radians)
%   I_top_of_t - cell-array [1 x n_mu] or [n_mu x 1] with arrays
%           for the electron fluxes at current energy as a function
%           of time. For the upward fluxes the natural condition is
%           to constrain the gradient of the upward fluxes to zero.
%           For example: {1e3*sin(10*2*pi*linspace(0,0.1,41)).^.2,...}
%           would be a suitable element for a harmonically flux
%           varying at 10 Hz from 0 to 0.1 s in 2.5 ms steps. All
%           elements must have matching sizes. For downward fluxes
%           the content should be all zero, e.g. {..., 0*ones(1,41)}.
%   I0    - initial electron fluxes [I0_mu1;I0_mu2;...;I0_mun]
%   v     - velocity (m/s), scalar double, typically v_of_E(E)
%   QC_mu - cascading contribution and local production
%           contributing to electron fluxes in each stream, double
%           array [(n_z*n_mu) x nt]
%   A     - total (electron collision cross section * neutral density),
%           double array [n_z x 1] (/m)
%   B_b2b - elastic collision from beams mu_j to beam mu_i, double
%           array [n_z x n_mu x n_mu] (/m)
%   Dz    - optional velocity diffusion parameter, defaults to
%           zero, slightly dubious use.
%   lowerBvals - depreceated input parameter, optional defaults to
%           1, unused regardless of value, was used for selecting
%           between setting fluxes or gradients of fluxes at lower
%           boundary to zero.
% Output:
%   Ie_zt - ( altitude x stream ) x time variation of electron
%           fluxes. 
% 
% Example: 
%   
% See also Ie_Mstream_tz_2_aurora Etrp_flickering_PnL Etrp_Flickering_PnLp8keV_9stream

%  Copyright © Bjorn Gustavsson 20180603, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later

% Input nargin-ruler (Handy to use tool for counting input arguments):
%                                       1    2  3          4  5 6     7 8     9 10 11         12


% The time-dependent Crank-Nicolson solution of the time-dependent
% 2-stream electron-transport equations 
% 
%  1/v*D/Dt(Id) = D/Dz(Id) - A*Id +B*Iu + QCd
%  
%  1/v*D/Dt(Iu) = -D/Dz(Iu) - A*Iu +B*Id + QCu
% 
% becomes: 
% 
%  Mlhs*Idu(t+dt) = Mrhs*Idu(t) + [bc(Id,z_min);QCd;bc(Id,z_max);bc(Iu,z_min);QCu;bc(Iu,z_max)]
% 
% Where Mr and Ml are built with band-diagonal and diagonal
% components.

if nargin < 11 || isempty(lowerBvals)
  lowerBvals = 1;
end

dZ = h_atm(2)-h_atm(1);
t = t_in;
dt = t(2)-t(1);
% Courant-Freidrich-Lewy number, should at least be small (<4)
C_CFL = v*dt/dZ;
n_factors = 2.^[0:22];
iFactor = 1;
while (1 < C_CFL) && (iFactor < 22)
  iFactor = iFactor+1;
  t = linspace(t_in(1),t_in(end),numel(t_in)*n_factors(iFactor)+1-n_factors(iFactor));
  dt = t(2)-t(1);
  C_CFL = v*dt/dZ;
end

C_CFLfactor = n_factors(iFactor);

%Ie_zt = zeros(numel(h_atm)*4,numel(t)); % *2 because [I_down;I_up]
if iscell(I0)
  Ie_zt = cell2mat(I0(:));
else
  Ie_zt = I0(:);
end

%Ad = diag(A + min(0,(ne_dLdE - neL*2/(sum(DE)))/DE(1)));
Ad = diag(A); % - min(0,(ne_dLdE - neL*2/(sum(DE)))/DE(1)));
% Bd = diag(B);

% Spatial differentiation diagonal matrix, for one stream
DdZ = diag(-1./(diff(h_atm)*4),-1) + diag(1./(diff(h_atm)*4),1);
% Temportal differentiation diagonal matrix, for one stream
Ddt = diag(1/(v*(t(2)-t(1)))*ones(size(h_atm)));


D2Zd = Dz/2*( diag(1./(diff(h_atm)).^2,-1) + ...
              diag(-2./(diff(h_atm([1,1:end]))).^2,0) + ...
              diag(1./(diff(h_atm)).^2,1) );

% Q_ee = repmat(neL,size(B_b2b,2),size(t,2)).*IeEpdE*2/(DE(2)*sum(DE));

%% 2-stream version Crank-Nicolson matrices
% These below work, (barring instabillity issues)
%
% Mrhs_d = [ DdZ+Ddt-Ad/2+D2Zd, Bd/2];
% Mrhs_u = [ Bd/2, -DdZ+Ddt-Ad/2+D2Zd];
%
% Mlhs_d = [-DdZ+Ddt+Ad/2-D2Zd, -Bd/2];
% Mlhs_u = [-Bd/2,  DdZ+Ddt+Ad/2-D2Zd];

Mlhs = [];
Mrhs = [];
for i2 = 1:size(B_b2b,2),
  MBlhs = [];
  MBrhs = [];
  for i3 = 1:size(B_b2b,2),
    if i2 ~= i3
      tmp_lhs = diag(-B_b2b(:,i2,i3)/2);
      tmp_lhs([1 end],:) = 0;
      tmp_rhs = diag(B_b2b(:,i2,i3)/2);
      tmp_rhs([1 end],:) = 0;
      
      MBlhs = [MBlhs,tmp_lhs]; % diag(-B_b2b(:,i2,i3)/2)];
      MBrhs = [MBrhs,tmp_rhs]; % diag(B_b2b(:,i2,i3)/2)];
    else
      tmp_lhs =  mu(i2)*DdZ+Ddt+Ad/2-D2Zd + diag(-B_b2b(:,i2,i3)/2);
      tmp_rhs = -mu(i2)*DdZ+Ddt-Ad/2+D2Zd + diag(B_b2b(:,i2,i3)/2);
      tmp_lhs(1,:) = 0; % fixed values in all beams at
      tmp_lhs(1,1) = 1; % lowest altitude
      if mu(i2) < 0 % downward fluxes, fixed values at top
        tmp_lhs(end,:) = 0;
        tmp_lhs(end,end) = 1;
      else % upward fluxes, zero gradients out at top
        tmp_lhs(end,:) = 0;
        tmp_lhs(end,end-1:end) = [-1,1];
      end
      tmp_rhs([1,end],:) = 0;
      MBlhs = [MBlhs,tmp_lhs];
      MBrhs = [MBrhs,tmp_rhs];
    end
  end
  Mlhs = [Mlhs;sparse(MBlhs)];
  Mrhs = [Mrhs;sparse(MBrhs)];
end

% F = factorize(A) ;  % computes the factorization of A
invMLHS = inverse(Mlhs);    % no flops, flags S as a factorized form of inv(A)

idx4QCbeams = sort([1:numel(h_atm):(numel(mu)*numel(h_atm)),...
    numel(h_atm):numel(h_atm):(numel(mu)*numel(h_atm))]);

I_zmax_of_t = cell2mat(I_top_of_t(:));
Ie_zt_prev = Ie_zt;
i_tin = 1;

for i_t = 2:(numel(t)),
  
  QC_I_zmax = ( QC_mu(:,i_tin)/2 + QC_mu(:,i_tin+1)/2 );
  Ibv = (I_zmax_of_t(:,i_tin)*[0,1])';
  QC_I_zmax(idx4QCbeams) = Ibv(:);
  
  % Crank-Nicolson step in time
  Ie_zt_next = invMLHS*([Mrhs]*Ie_zt_prev + QC_I_zmax);
  Ie_zt_prev = Ie_zt_next;
  
  % if C_CFLfactor == 1 || rem(i_t,n_factors(iFactor))==1
  if C_CFLfactor == 1 || rem(i_t,C_CFLfactor)==1
    i_tin = i_tin+1;
    Ie_zt(:,i_tin) = Ie_zt_next;
  end

end
% $$$ 
% $$$ for i_t = 2:(numel(t)),
% $$$   
% $$$   % QC_I_zmax are built with the pattern:
% $$$   % [I_e_mu1_bc(zmin);QC_mu1(it+1/2);Ie_mu1_bc(zmax);...]
% $$$   % pitch-angle-stream by pitch-angle-stream...
% $$$ % $$$   QC_I_zmax = ( QC_mu(:,i_t)/2 + QC_mu(:,i_t+1)/2 + ...
% $$$ % $$$                 Q_ee(:,i_t)/2  + Q_ee(:,i_t+1)/2 );
% $$$   QC_I_zmax = ( QC_mu(:,i_t)/2 + QC_mu(:,i_t+1)/2 );
% $$$   Ibv = (I_zmax_of_t(:,i_t)*[0,1])';
% $$$   QC_I_zmax(idx4QCbeams) = Ibv(:);
% $$$ 
% $$$   % Crank-Nicolson step in time
% $$$   Ie_zt(:,i_t+1) = invMLHS*([Mrhs]*Ie_zt(:,i_t) + QC_I_zmax);
% $$$   
% $$$ end
