function [LHS,RHS] = de_M_stream_checker(Ie_zt,h_atm,t_in,mu,I_top_of_t,I0,v,QC_mu,A,B_b2b,Dz,lowerBvals)
% de_M_stream_CNztusD - multistream Crank-Nicholson electron transport
% solver with up-stream differential operator for spatial
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
%  [Ie_zt] = de_M_stream_CNztus(h_atm,t,mu,I_top_of_t,I0,v,...
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
%                                         1    2  3          4  5 6     7 8     9 10 11         12


% The time-dependent 2-stream electron-transport equations 
% 
%  D/Dt(Id) = -mu_avgd*D/Dz(Id) - A*Id +B*Iu + QCd
%  
%  D/Dt(Iu) = -mu_avgu*D/Dz(Iu) - A*Iu +B*Id + QCu
% 
% becomes: 
% 
%  1/v*[dId_dt;dIu_dt] + [mu_avgd*dId_dz;mu_avgu*dIu_dz] = [A+B11,B21;B12,A+B22]*[Id;Iu] + [QCd;QCu]
% 

if nargin < 11 || isempty(lowerBvals)
  lowerBvals = 1;
end

dZ = h_atm(2)-h_atm(1);

t = t_in;

C_CFLfactor = n_factors(max(1,iFactor-1));
%disp(n_factors(iFactor))
if iscell(I0)
  Ie_zt = cell2mat(I0(:));
else
  Ie_zt = I0(:);
end

Ad = diag(A)/2 + diag(A([2:end,end]))/2; % Shift half-step up
Au = diag(A)/2 + diag(A([1,1:end-1]))/2; % shift half-step down
% Bd = diag(B);

% Spatial differentiation diagonal matrix, for one stream
% Here we tuck on a fictious height one tick below lowest height
% just to make it possible to use the diff function.
h4diff = [h_atm(1) - (h_atm(2)-h_atm(1));h_atm];

% Up-stream differential operators
DdZd = diag(-1./(diff(h4diff)*2),0) + diag(1./(diff(h4diff(1:end-1))*2),1);
DdZu = diag(-1./(diff(h4diff(1:end-1))*2),-1) + diag(1./(diff(h4diff)*2),0);

% Temporal differentiation diagonal matrix, for one stream
Ddt = diag(1/(v*(t(2)-t(1)))*ones(size(h_atm)));

% Diffusion operator, typically all zeroes.
D2Zd = d2M(h_atm);

D2Zd(1,1) = 0; % Will allways be OK since first row will be set to
               % zeroes with one one at the right diagonals to keep
               % fixed lower-boundary-values.

%% 2-stream version Crank-Nicolson matrices
% The matrices for the Crank-Nicolson scheme are simple enough for
% a 2-stream system with only losses (due to inelastic collisions
% degrading electrons out of beam and elastic back-scattering to
% the other stream). These below work. 
%
% Mrhs_d = [ DdZ+Ddt-Ad/2+D2Zd,          Bd/2     ];
% Mrhs_u = [         Bd/2,      -DdZ+Ddt-Au/2+D2Zd];
%
% Mlhs_d = [-DdZ+Ddt+Ad/2-D2Zd,         Bd/2     ];
% Mlhs_u = [        -Bd/2,      DdZ+Ddt+Au/2-D2Zd];

% For the general case we have a larger number of streams with
% elastic scattering between them, this makes the system larger
% with a larger number of blocks with diagonals for the elastic
% stream-to-stream scattering, but the pattern remains.
Mlhs = [];
Mrhs = [];
for i2 = 1:size(B_b2b,2),
  MBrhs = [];
  for i3 = 1:size(B_b2b,2),
    if mu(i2) < 0 % down-ward fluxes, shift half-step up
      Bb2b = B_b2b(:,i2,i3)/2 + B_b2b([2:end,end],i2,i3)/2;
    else          % up-ward fluxes, shift half-step down
      Bb2b = B_b2b(:,i2,i3)/2 + B_b2b([1,1:end-1],i2,i3)/2;
    end
    if i2 ~= i3
      %tmp_lhs = diag(-B_b2b(:,i2,i3)/2);
      tmp_rhs = diag(Bb2b);
      tmp_rhs([1 end],:) = 0;
      
      MBrhs = [MBrhs,tmp_rhs];
    else
      C_D = Dz(min(i2,numel(Dz)));
      if mu(i2) < 0 % down-ward fluxes,
        DdZ = DdZd;
        Ac = Ad;
      else % up-ward fluxes
        DdZ = DdZu;
        Ac = Au;
      end
      tmp_rhs = -Ac+C_D*D2Zd + diag( Bb2b);
      tmp_rhs([1,end],:) = 0;
      MBrhs = [MBrhs,tmp_rhs];
    end
  end
  Mrhs = [Mrhs;sparse(MBrhs)];
end

nZ = numel(h_atm);
%% Somewhere(Here?) calculate d(Ie_ztE)/dt and d(Ie_ztE)dz
for i_mu = numel(mu):-1:1,
  [dIedt((1:nZ)+nZ*(i_mu-1),:),mu_dIedz((1:nZ)+nZ*(i_mu-1),:)] = gradientC(Ie_ztE(:,:,iE),t,h_atm);
  mu_dIedz((1:nZ)+nZ*(i_mu-1),:) = mu(i_mu)*mu_dIedz((1:nZ)+nZ*(i_mu-1),:);
end
LHS = 1/v*dIedt + mu_dIedz;



QC_I_zmax = ( QC_mu(:,[1:end])/2 + QC_mu(:,[2:end,end])/2 );
QC_I_zmax(idx4QCbeams) = Ibv(:);

% Crank-Nicolson step in time
RHS =[Mrhs]*Ie_zt + QC_I_zmax;
