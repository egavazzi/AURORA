function [Ie_z] = de_M_stream_us(h_atm,mu,I_top_of_mu,QC_mu,A,B_b2b,lowerBvals)
% de_M_stream_CNztus - multistream electron transport solver
% solver with up-stream differential operator for spatial
% derivatives. This integrates the steady-state electron
% transport equation:
% 
%   mu*DI_mu/Dz = -A*I_mu + B*I_mu' + Q+C
% 
% For electron fluxes at energy E in an arbitrary number of
% pitch-angle streams with average cosine-of-pitch-angle mu. The
% upper boundary condition I_top_of_mu is to be a cell-array with
% the upper boundary condition as a function of time pitch-angle by
% pitch-angle.
% 
% Calling:
%  [Ie_z] = de_M_stream_CNztus(h_atm,t,mu,I_top_of_mu,I0,v,...
%                               QC_mu,A,B,Pmu2mup,Dz,...
%                               lowerBvals)
% Input:
%   h_atm - distance along themagnetic field (tilted altitude) [nz x 1] (m)
%   mu    - pitch-angle [1 x n_mu] (radians)
%   I_top_of_mu - cell-array [1 x n_mu] or [n_mu x 1] with arrays
%           for the electron fluxes at current energy as a function
%           of time. For the upward fluxes the natural condition is
%           to constrain the gradient of the upward fluxes to zero.
%           For example: {1e3*sin(10*2*pi*linspace(0,0.1,41)).^.2,...}
%           would be a suitable element for a harmonically flux
%           varying at 10 Hz from 0 to 0.1 s in 2.5 ms steps. All
%           elements must have matching sizes. For downward fluxes
%           the content should be all zero, e.g. {..., 0*ones(1,41)}.
%   QC_mu - cascading contribution and local production
%           contributing to electron fluxes in each stream, double
%           array [(n_z*n_mu) x nt]
%   A     - total (electron collision cross section * neutral density),
%           double array [n_z x 1] (/m)
%   B_b2b - elastic collision from beams mu_j to beam mu_i, double
%           array [n_z x n_mu x n_mu] (/m)
%   lowerBvals - depreceated input parameter, optional defaults to
%           1, unused regardless of value, was used for selecting
%           between setting fluxes or gradients of fluxes at lower
%           boundary to zero.
% Output:
%   Ie_z - ( altitude x stream ) x time variation of electron
%           fluxes. 
% 
% Example: 
%   
% See also Ie_Mstream_tz_2_aurora Etrp_flickering_PnL Etrp_Flickering_PnLp8keV_9stream

%  Copyright © Bjorn Gustavsson 20180603, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later

% Input nargin-ruler (Handy to use tool for counting input arguments):
%                                         1    2  3          4  5 6     7 8     9 10 11         12


% Pattern of 2-stream transport matrix is:
%
%   [ A+D , -B  ]   [I_d]
%   [           ] * [   ]
%   [ -B  , A+D ]   [I_u]
% 
%[                                                                               ] [Idhm1]
%[             0, 1/dH + Ah/2,-1/dH + Ahp1/2,           0,     -Bh/2,     -Bhp1/2] [  Idh]
%[                                                                               ] [Idhp1]
%[                                                                               ]*[Iuhm1]
%[       -Bhm1/2,       -Bh/2,             0, -1/dH + Ahm1, Ah + 1/dH,          0] [  Iuh]
%[                                                                               ] [Iuhp1]
%
%
%
%[ Ahm1/2 + 1/dH, Ah/2 - 1/dH,             0,     -Bhm1/2,     -Bh/2,           0] [Idhm1]
%[             0, Ah/2 + 1/dH, Ahp1/2 - 1/dH,           0,     -Bh/2,     -Bhp1/2] [  Idh]
%[             0,           0,             1,           0,         0,           0] [Idhp1]
%[             0,           0,             0,           1,         0,           0]*[Iuhm1]
%[       -Bhm1/2,       -Bh/2,             0, Ahm1 - 1/dH, Ah + 1/dH,           0] [  Iuh]
%[             0,       -Bh/2,       -Bhp1/2,           0, Ah - 1/dH, Ahp1 + 1/dH] [Iuhp1]
%
%
%

if nargin < 11 || isempty(lowerBvals)
  lowerBvals = 1;
end
% Ad = diag(A);
Ad = diag(A/2+A([2:end,end])/2); % Shifted half-step up
Au = diag(A/2+A([1,1:end-1])/2); % Shifted half-step down  

dZ = h_atm(2)-h_atm(1);

h4diff = [h_atm(1) - (h_atm(2)-h_atm(1));h_atm];

% Up-stream differential operators
DdZd = diag(-1./(diff(h4diff)*2),0) + diag(1./(diff(h4diff(1:end-1))*2),1);
DdZu = diag(-1./(diff(h4diff(1:end-1))*2),-1) + diag(1./(diff(h4diff)*2),0);

Mlhs = [];
for i2 = 1:size(B_b2b,2),
  MBlhs = [];
  for i3 = 1:size(B_b2b,2),
    if i2 ~= i3
      if mu(i2) < 0 % down-ward fluxes,
        % Shift half-step up
        Bb2b = diag(-(B_b2b(:,i2,i3)/2 + B_b2b([2:end,end],i2,i3)/2));
      else % up-ward fluxes
        % shift half-step down
        Bb2b = diag(-(B_b2b([1,1:end-1],i2,i3)/2 + B_b2b(:,i2,i3)/2));% diag(-B_b2b(:,i2,i3));
      end
      tmp_lhs = Bb2b; % diag(-B_b2b(:,i2,i3));
      tmp_lhs([1 end],:) = 0;
      MBlhs = [MBlhs,tmp_lhs];
    else
      if mu(i2) < 0 % down-ward fluxes,
        % A and B_b2b should be shifted half a step up
        DdZ = DdZd;
        Ac = Ad;
        Bb2b = diag(-(B_b2b(:,i2,i3)/2 + B_b2b([2:end,end],i2,i3)/2));
     else % up-ward fluxes
        % A and B_b2b should be shifted half a step down
        DdZ = DdZu;
        Ac = Au;
        Bb2b = diag(-(B_b2b([1,1:end-1],i2,i3)/2 + B_b2b(:,i2,i3)/2));% diag(-B_b2b(:,i2,i3));
      end
      tmp_lhs =  mu(i2)*DdZ+Ac + Bb2b; % diag(-B_b2b(:,i2,i3));
      tmp_lhs(1,:) = 0; % fixed values in all beams at
      tmp_lhs(1,1) = 1; % lowest altitude
      if mu(i2) < 0 % downward fluxes, fixed values at top
        tmp_lhs(end,:) = 0;
        tmp_lhs(end,end) = 1;
      else % upward fluxes, zero gradients out at top
        tmp_lhs(end,:) = 0;
        tmp_lhs(end,end-1:end) = [-1,1];
      end
      MBlhs = [MBlhs,tmp_lhs];
    end
  end
  Mlhs = [Mlhs;sparse(MBlhs)];
end

% F = factorize(A) ;     % computes the factorization of A
invMLHS = inverse(Mlhs); % no flops, flags S as a factorized form of inv(A)

idx4QCbeams = sort([1:numel(h_atm):(numel(mu)*numel(h_atm)),...
    numel(h_atm):numel(h_atm):(numel(mu)*numel(h_atm))]);

I_zmax_of_t = cell2mat(I_top_of_mu(:));
% Ie_z_prev = Ie_z;
% i_tin = 1;


QC_I_zmax = ( QC_mu );
Ibv = (I_zmax_of_t*[0,1])';
QC_I_zmax(idx4QCbeams) = Ibv(:);

Ie_z = invMLHS*QC_I_zmax;
