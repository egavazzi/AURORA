function Brratio = branching_resample(WL,wlbrratio)
% BRANCHING_RESAMPLE resapmple the photo-ionization branching ratios. 
%
% Calling:
%  Brratio = branching_resample(WL,wlbrratio)
% Input:
%  WL - wavelenths array to resample to. 
%  WLBRATIO - table of threshold wavelengths and branching ratios.
% Output:
%  Brration - branching ration between states of ions after
%  photo-ionization

% Copyright: B. Gustavsson 20100527

Brratio = zeros([length(WL) size(wlbrratio,2)]);
for i = size(wlbrratio,1):-1:1,
  
  for j = 2:size(wlbrratio,2),                     
    
    Brratio(WL<wlbrratio(i,1),j) = wlbrratio(i,j);
    
  end
  
end
Brratio(:,1) = WL';
