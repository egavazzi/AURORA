%% Animate Auroral Rays
% This is a script utilizing the AIDA-tools toolbox for 
% projecting 3-D volume emission-rates down to an image-
% plane. The script is unfortunately not self-contained,
% Please contact the author for additional information.

%   Copyright © 2018-2019 Bjorn Gustavsson, <bjorn.gustavsson@uit.no>
%   This is free software, licensed under GNU GPL version 2 or later

root_dir = '/media/bgu001/5f5e8978-a828-4fd4-aabf-2032a3fb895b/Data/Bjorns-panic-repository/Etrp-tdms/Run20190619/';
load_dirs = {'Awa_1200-1'
             'Awa_1200-2'
             'Awa_1200-4'
             'AWA_2400-1'
             'AWA_2400-2'
             'AWA_2400-4'
             'AWA_3000-2'
             'AWA_3000-4'};


for iD = 1:numel(load_dirs)
  cd(fullfile(root_dir,load_dirs{iD}))
  pwd
  load Qzt_all_L t h_atm Q4278 Q6730 Q7774 Q8446 QO1S QO1D
  I5577l = 0*Q4278;
  for it = 2:numel(t),
    I5577l(:,it) = I5577l(:,it-1)*exp(-(t(it)-t(it-1))/0.7) + QO1S(:,it);
  end
  Ibeams(iD).t = t;                                  
  Ibeams(iD).z = h_atm/1e3;
  Ibeams(iD).wavelengths = [4278 5577 6730 7774 8446];
  Ibeams(iD).I{1} = Q4278;                              
  Ibeams(iD).I{2} = I5577l;                             
  Ibeams(iD).I{3} = Q6730;
  Ibeams(iD).I{4} = Q7774;
  Ibeams(iD).I{5} = Q8446;
  
end



idx_col_ray = [30 42 54 66];
idx_row_ray = [25 75];
[idx_col_ray,idx_row_ray] = meshgrid(idx_col_ray,idx_row_ray);
idxCR = [idx_col_ray(:),idx_row_ray(:)]

I3D = 0*ZfI;
z3D = squeeze(ZfI(1,1,:));

zi = squeeze(ZfI(1,1,:));
wbh = waitbar(1,'patience Bjeorn, patience...');
for it = 151:-1:1,
  
  I3D = 0*ZfI;
  for i_L = 1:numel(Ibeams(1).I),
    for iRay = 1:numel(Ibeams)
      t = Ibeams(iRay).t;
      z = Ibeams(iRay).z;
      % if it <= numel(t)
        I3D(idxCR(iRay,2),idxCR(iRay,1),:) = ...
            interp1(z,Ibeams(iRay).I{i_L}(:,min(it,end)),zi,'linear',0);
        I_of_z(:,iRay) = Ibeams(iRay).I{i_L}(:,min(it,end));
      % end
    end
    if exist('fKray')
      I3D = convn(I3D,fKray,'same');
    end
    Imgstacks{i_L}(:,:,it) = fastprojection(I3D,...
                                            stns.uv,stns.d,stns.l_cl,stns.bfk,256*[1 1],...
                                            ones(256,256),stns.sz3d);
    
  end
  try
    waitbar(it/151,wbh)
    if it == 50
      waitbar(it/151,wbh,'2 3rds to the end now')
    elseif it == 49
      waitbar(it/151,wbh,'soon it will be done...')      
    end
  end
end

try
  close(wbh)
end