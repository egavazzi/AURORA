function [varargout] = get_all_photo_ion_xs(name,lambda)
% GET_ALL_PHOTO_ION_XS - collect all photo-ionization cross sections
% All cross sections for species 'NAME' at wavelengths LAMBDA (in m)
% NAME should be the name of a species 'NO', 'N2' or 'O2',
% for which to calculate the electron cross section. LAMBDA -
% wavelength of photons in m. 
% 
% Calling:
%  [XS,E_threshold,level_name,Brancing_ratio] = get_all_photo_ion_xs(name,lambda);
% Input:
%  name - name of species,{'NO','O2','N2'}
%  lamba    - Photon wavelength vector [1 x n] (m)
% Output:
%  Xs - photo-ionization cross sections (m^2). The
%       photo-ionization cross sections are sorted%

% Copyright Bjorn Gustavsson 20100527

c0	= 2.99792458e8;		    % Speed of light [m/s]
h	= 6.62618e-34;		    % Plank's constant [Js]
q_e     = 1.6021773e-19;            % elementary charge [C]

fid = fopen([name,'_ion_levels.dat'],'r');

C = textscan(fid,'%f %s %s','CommentStyle','%');
fclose(fid);

E_threshold = C{1};

E = (h*c0/q_e)./lambda;
XS = [];
iXS = 1;
for iLevel = 1:length(E_threshold),
  
  % disp((['p_',C{2}{iLevel}]))
  %try
    load(['p_',C{2}{iLevel}])
    xs1 = interp1(EnXS(:,1),EnXS(:,2),E,'pchip');
    xs2 = min(EnXS(end,2),exp(interp1(EnXS(:,1),log(EnXS(:,2)),E,'linear','extrap')));
    XS(iXS,:) = xs1;
    XS(iXS,E>max(EnXS(:,1))) = xs2(E>max(EnXS(:,1)));
    XS(iXS,E<E_threshold(iLevel)) = 0;
    XS = max(0,XS);
    E_Threshold(iXS) = E_threshold(iLevel);
    Level{iXS} = C{2}{iLevel};
    % disp(['Loading photo-ionization cross section for: ',Level{iXS}])
    % plot(E,XS(iXS,:),'b')
    % hold on
    % plot(EnXS(:,1),EnXS(:,2),'r.')
    % title(C{2}{iLevel})
    % pause(0.1)
    % hold off
    iXS = iXS + 1;
  %catch
  %  disp(['Failing to make photo-ionization cross section for: ',C{2}{iLevel}])
  %end
  
end

XS(~isfinite(XS(:))) = 0;
if nargout > 0
  % Ionization Cross sections
  varargout{1} = XS;
end
if nargout > 1
  % Energy levels of ion states
  varargout{2} = E_Threshold;
end
if nargout > 2
  % Ion state names/labels
  varargout{3} = Level;
end
if nargout > 3
  % Branching ratios
  Branching_ratio = XS(2:end,:)./repmat(sum(XS(2:end,:)),size(XS(2:end,:),1),1);
  Branching_ratio(~isfinite(Branching_ratio(:))) = 0;
  varargout{4} = Branching_ratio;
  
end

