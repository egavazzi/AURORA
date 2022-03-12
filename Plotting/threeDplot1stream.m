%% Script producing a 3-cuts plot
% of differential electron-fluxes in one (the most field-aligned)
% electron-stream.
% 
% Before runnign this script electron-fluxes needs to be read
% in. This is typically done with:
% 
%  [t,h_atm,E,mu_lims,IeZTE,mu_scatterings] = Ie_ztE_loader({'PnLs-07-4000-3keV-0-10'});
%  dE = diff(E);dE(end+1) = dE(end);

%% t-E-z-cube-corners
[rCtEx,rCtEy,rCtEz] = meshgrid(t([1 end]),E([1 size(IeZTE,3)]),h_atm([1 end])/1e3); 
rC = [rCtEx(:),rCtEy(:),rCtEz(:)];

%% Slice plot
slice(t,E(1:size(IeZTE,3)),h_atm/1e3,...
      log10(permute(IeZTE(1:numel(h_atm),:,:),[3 2 1]) ./ ...
            repmat(dE(1:size(IeZTE,3))',[1,numel(t),numel(h_atm)])),...
      0.5,2500,150),
shading flat
axis tight
caxis(max(caxis)+[-5 0])
set(gca,'ydir','normal','yscale','log')                        
set(gca,'ytick',[10 100 1000],'YMinorGrid','off','fontsize',12)
xlabel('time (s)','fontsize',14)
ylabel('Energy (eV)','fontsize',14)
zlabel('height (km)','fontsize',14)
cblh = colorbar_labeled('#e^-/m^2/s/eV','log');
set(cblh,'position',get(cblh,'position')+[-0.01 0.1 0 -0.2])
set(cblh,'position',get(cblh,'position')+[-0.0 0.1 0 -0.2]) 
set(cblh,'position',get(cblh,'position')-[-0.0 0.1 0 -0.2])
set(cblh,'position',get(cblh,'position')+[-0.0 0.05 0 -0.1])
set(cblh,'position',get(cblh,'position')+[-0.0 0.01 0 -0.1])
set(cblh,'position',get(cblh,'position')+[-0.0 0.0 0 0.05])  
set(cblh,'position',get(cblh,'position')-[-0.0 0.0 0 0.01])
hold on
plotCube(rC([1 3 4 2 5 7 8 6],:),'linewidth',2,'color','k')                
plot3(rC([1 3 4 2 1],1),...
      rC([1 3 4 2 1],2),...
      150*ones(5,1),...
      'g','linewidth',2) 
plot3(0.5*ones(5,1),...
      rC([1 5 6 2 1],2),...
      rC([1 5 6 2 1],3),...
      'color',[0 1 1],'linewidth',2) 
plot3(rC([1 3 7 5 1],1),...
      2500*ones(5,1),...
      rC([1 3 7 5 1],3),...
      'm','linewidth',2)
set(cblh,'position',get(cblh,'position')-[-0.05 0.0 0 0.0])
set(cblh,'position',get(cblh,'position')-[0.05 0.0 0 0.0]) 
set(cblh,'position',get(cblh,'position')-[0.02 0.0 0 0.0])
title('Electron flux in stream #1, \theta: 0-10^\circ','fontsize',14)
set(cblh,'position',get(cblh,'position')-[0.0 0.01 0 0.0]) 
set(cblh,'position',get(cblh,'position')-[0.0 -0.005 0 0.0])
set(cblh,'position',get(cblh,'position')-[0.0 0.0025 0 0.0]) 

plot3(rC([1 3],1),2500*ones(2,1),[150 150],...
      'color',[1 1 1]/2,'linewidth',2)
plot3([0.5 0.5],2500*ones(2,1),rC([1 5],3),...
      'color',[1/2 1/2 1],'linewidth',2)
plot3([0.5 0.5],rC([1 2],2),[150 150],...
      'color',[0 1 1/2],'linewidth',2)


%   Copyright © 2019 Bjorn Gustavsson, <bjorn.gustavsson@irf.se>
%   This is free software, licensed under GNU GPL version 2 or later
