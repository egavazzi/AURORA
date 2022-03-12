function [varargout] = photo_e_source(I_p,wl,Eout,varargin)
% PHOTO_E_SOURCE - Source of photo-ionized electrons and ions
% 
% Calling:
%   pho_e_flux = photo_e_source(I_p,wl,Eout,n_j,s_i_j,b_j,...)
%   [pho_e_flux,ion_rate_j,...] = photo_e_source(I_p,wl,Eout,n_j,s_i_j,b_j,...)
%   pho_e_flux = photo_e_source(I_p,wl,Eout,'demo',n_j,s_i_j,b_j,...)
%   [pho_e_flux,ion_rate_j,...] = photo_e_source(I_p,wl,Eout,'demo',n_j,s_i_j,b_j,...)
% 
% Input: 
%   I_p  - solar photon flux [n x m] array of photons (m^{-2}s^{-1}). 
%   WL   - Wavelength [1 x m] array defining the wavelength bins.
%   Eout - the requested energies of the output [1 x M] (eV).
%   'demo' - optional flag, if used the function will display the
%            production of secondary electrons wavelength for
%            wavelength, and species by species.
%   N_j  - number density [n x 1] of atmospheric species j (m^{-3}).
%   S_j  - ionisation cross section (m^{-2}).
%   B_j  - is the branching ratio matrix with wl in column 1 and
%          branchin factors in the following. 
%   The function handles arbitrary number of species provided that
%   density, ionization cross section and branching ratios are
%   given for all species.
% Output:
%   pho_e_flux - photo-electron production rate [n x M] (#/m^3/s/eV?)
%
% NOTE!
%  There is some uncertainty whether this function expects photon
%  flux in (/A) or (/gradient(WL)), the same uncertainty is for the
%  output (/eV) or (/gradient(Eout))

% Copyright: B. Gustavsson 20100527

% photo_e_source - calculates the source of photoelectrons.
% wavelengths, N1, N2, N3 are 

% Setting up physical constants
c0	= 2.99792458e8;		    % Speed of light [m/s]
h_planck= 6.62618e-34;		    % Plank's constant [Js]
kB	= 1.380662e-23;		    % Boltzmann constant [J/K]
q_e     = 1.6021773e-19;            % elementary charge [C]

sI = size(I_p);
nwl = sI(2);
nh = sI(1);

pho_e_flux = zeros([nh,length(Eout)]);

if strcmp(varargin{1},'demo')
  
  varargin = varargin(2:end);
  demo = 1;
  
else
  
  demo = 0;
  
end

if length(varargin)
  
  for i1=1:3:length(varargin),
    
    n = varargin{i1};
    s_i = varargin{i1+1};
    branch = varargin{i1+2};
    %ionization_rate_i = 0;
    
    for i2 = 1:(size(branch,2)-1)
      
      % [/m^3/s]         [/m^3]        [#/m^2/s]      [m^2]
      ion_rate = ( repmat(n,[1 nwl]) .* I_p .* repmat(s_i'.*branch(:,i2+1)',[nh 1]) );
      % This is, somewhat surprisingly, faster by about a factor of
      % 2 (matlab 2013a) than the variant below with an outer
      % product between density profile and cross-section:
      % ion_rate = I_p .* (n* (s_i'.*branch(:,i2+1)'));
      
      i_0plus = find(branch(:,i2+1)>0);
      i_last = i_0plus(end);
      wl_th = branch(i_last,1);
      
      I = find(~isfinite(ion_rate(:))|ion_rate(:)<0);
      ion_rate(I) = 0;
      
      %ionization_rate_i = ionization_rate_i + ion_rate;
      ionization_rate_i(:,i2) = sum(ion_rate,2);
      
      E = h_planck*c0/q_e*(1./wl-1/wl_th);
      %pho_e_flux = pho_e_flux + interp1(E,ion_rate',Eout,'pchip')';
      
      pho_e_flux(:,1) = pho_e_flux(:,1) + sum(ion_rate(:,E<(Eout(1)+Eout(2))/2),2);
      for i3 = 2:length(Eout)-1
        pho_e_flux(:,i3) = pho_e_flux(:,i3) + sum(ion_rate(:,(Eout(i3-1)+Eout(i3))/2<E&E<(Eout(i3)+Eout(i3+1))/2),2);
      end
      pho_e_flux(:,end) = pho_e_flux(:,i3) + sum(ion_rate(:,E>(Eout(i3)+Eout(i3+1))/2),2);
      
      if demo
        
        pcolor(Eout,1:size(I_p,1),log10(pho_e_flux))
        shading flat,colorbar,axis xy,cax = caxis;caxis([-6 0]+cax(2))
        title('log_{10} photoelectron source spectra (m^{-3}s{-1}eV^{-1})','fontsize',16)
        xlabel('electron energy (eV)','fontsize',16)
        drawnow
      end
      
    end
    %disp([length(varargin) nargout i1 (i1-1)/3+2])
    if (i1-1)/3+2 <= nargout
      varargout{(i1-1)/3+2} = ionization_rate_i;
    end
  end
  
end

varargout{1} = pho_e_flux;
