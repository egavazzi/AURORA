function [OptDepth] = opt_depth(s_abs,phi,h,n,hout)
% OPT_DEPTH is a function that calculates the optical depth
% at the altitudes H and solar zenith angle PHI given the
% atmospherical spieces model N and the absorbtion cross sections
% S_ABS.
%
% Calling
% [OptDepth] = opt_depth(s_abs,phi,h,n,hout)
%
% Input:
%  S_ABS - Absorption cross section (m^{-2}) [NxM] array with
%          M being the number of species in the atmosphere and M
%          the number of wavelength intervalls
%  PHI   - solar zenith angle (radians) [1x1]
%  H     - altitude (m) array, [Lx1]
%  N     - Atmospheric concentrations (m^{-3}), [LxM]
%  hout  - {Optional} output altitude (m) [1x1]
%
% Output
%  OptDepth optical depth
%
% After M. H. Rees 1989

% Copyright B. Gustavsson 2006 0502

%re = 6.371e3;
re    = 6.378e6;            % Radius of earth [m]

if numel(phi)>1
  
  error('Size of PHI way too large - only work with 1 solar zenith angle at a time')
  
end
if nargin==5
  
  if numel(hout)>1
    
    error('Size of HOUT too large...')
    
  end
  
end

if phi <= pi/2
  
  phi = min(phi,pi/2-5e-3);
  OptDepth = zeros(length(h),length(s_abs(:,1)));
  for i_z = 1:(length(h)-1),
    
    for i_n = 1:size(n,2)
      
      n_i = n(:,i_n);
      maxih = length(h);
      
% $$$       integrand = n_i(i_z:end) .* (1-((re+h(i_z))./(re+h(i_z:end))).^2*sin(phi)^2).^(-.5).*gradient(h(i_z:end));
% $$$       cd2 = sum(integrand);
      integrand = n_i(i_z:end) .* (1-((re+h(i_z))./(re+h(i_z:end))).^2*sin(phi)^2).^(-.5);
      cd2 = trapz(h(i_z:end),integrand);
      
      part_depth = cd2' * s_abs(:,i_n)';
      OptDepth(i_z,:) = OptDepth(i_z,:) + part_depth;
      
    end
    
  end
  
else
  % When solar zenith angle > 90 degrees geometry becomes different!  
  if nargin > 4 && ~isempty(hout) % exist('hout')
    H = sort([h;hout]);
  else
    H = sort(h);
  end
  
  z = min(H);
  minh = (re+z)*sin(pi-phi)-re;
  h_pi_2 = (re+h)*sin(pi-phi)-re;
  for i_n = 1:min(size(n))
    N_hor(:,i_n) = exp(interp1(h,log(n(:,i_n)),h_pi_2,'linear','extrap'));
  end
  
  ODhorizon = opt_depth(s_abs,pi/2,h_pi_2,N_hor);
  OptDepth  = opt_depth(s_abs,pi-phi,h,n);
  
% $$$   for i_n = 1:min(size(n))
% $$$     
% $$$     n_i = n(:,i_n);
% $$$     H = [minh:diff(h(1:2)):h(1)];
% $$$     
% $$$     H = unique([H';h]);
% $$$     N_i = exp(interp1(h,log(n_i),H,'linear','extrap'));
% $$$     
% $$$     maxih = length(H);
% $$$     
% $$$     integrand = N_i(1:maxih) .* (1-((re+H(1))./(re+H(1:maxih))).^2*sin(phi)^2).^(-.5).*gradient(H(1:maxih));
% $$$     % integrand(h_pi_2<70e3) = integrand(h_pi_2<70e3)*1e6;
% $$$     cd2 = flipud(cumsum(flipud(integrand)))';
% $$$     
% $$$     part_depth = cd2' * s_abs(:,i_n)';
% $$$     
% $$$     if ( i_n == 1 )
% $$$       
% $$$       OptDepth = part_depth;
% $$$       
% $$$     else
% $$$       
% $$$       OptDepth = OptDepth + part_depth;
% $$$       
% $$$     end
% $$$     
% $$$   end
% $$$   OptDepth = interp2(1:length(s_abs),H,OptDepth,1:length(s_abs),h);
  OptDepth = 2*ODhorizon - OptDepth;
  
end
