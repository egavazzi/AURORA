%% Reproduce Differential-flux-plots of Peticolas and Lummerzheim
% This script reproduces the differential electron flux panels of
% Figures 5 and 6 in Peticolas and Lummerzheim (2000), with fluxes
% from AURORA. This script requires the results from a 0.7 s long
% integration of 5 Hz square-wave modulated FAB from 4000 km. This
% is generated with one of the parameter-sets in
% Loop_finalRun_PnL_flickering10beams.m, that will be automatically
% saved in a directory named PnLs-05-4000-3keV-0-10


try
  [t_PnL,h_atm,E,mu_lims,IeZTE] = Ie_ztE_loader({'PnLs-05-4000-3keV-0-10'});
catch
  [t_PnL,h_atm,E,mu_lims,IeZTE] = Ie_ztE_loader({'PnLs-05-4000-3keV-FA'});
end

szIzte = size(IeZTE);
dE = diff(E);
dE(end+1) = dE(end);
[dZ490,iz490] = min(abs(490-h_atm/1e3))
[dz153,iz153] = min(abs(153-h_atm/1e3))

Ie_ZTE490 = squeeze(sum(IeZTE(iz490:numel(h_atm):end,:,:)));
Ie_ZTE153 = squeeze(sum(IeZTE(iz153:numel(h_atm):end,:,:)));

figure
sp211 = subplot(2,1,1);
pcolor(t_PnL,E(1:szIzte(3)),log10(Ie_ZTE490.'./repmat(dE(1:szIzte(3))'*4*pi,size(t_PnL)))),shading flat
caxis([-5 0]+12-1)
set(gca,'fontsize',12,'yscale','log','tickdir','out','box','off')
title('Electron flux, z = 490 km','fontsize',16)
ylabel('Energy (eV)','fontsize',15)
% cblh = colorbar_labeled('Log_{10}(Flux(m^{-2}eV^{-1}s^{-1}))','log','fontsize',14);
cblh1 = colorbar_labeled('Flux m^{-2}eV^{-1}s^{-1})','log','fontsize',14);
set(cblh1,'position',get(cblh1,'position')+[-0.01 0 0 0])
set(cblh1,'position',get(cblh1,'position')+[-0.01 0 0 0])
set(cblh1,'position',get(cblh1,'position')+[-0.01 0 0 0])
set(cblh1,'position',get(cblh1,'position')+[-0.01 0 0 0])
set(cblh1,'position',get(cblh1,'position')+[-0.01 0 0 0])
set(gca,'position',get(gca,'position')+[-0.01 0 0 0])
set(gca,'position',get(gca,'position')+[-0.01 0 0 0])
set(cblh1,'ytick',[1e6 1e8 1e10 1e12])
set(gca,'xticklabel','')
sp211pos = get(gca,'position');
set(gca,'ytick',[10 100 1000])

sp212 = subplot(2,1,2);
pcolor(t_PnL,E(1:szIzte(3)),log10(Ie_ZTE153.'./repmat(4*pi*dE(1:szIzte(3))',size(t_PnL)))),shading flat
set(gca,'fontsize',12,'yscale','log','tickdir','out','box','off')
caxis([-5 0]+12)
title('Electron flux, z = 153 km','fontsize',16)
xlabel('time (s)','fontsize',15)
ylabel('Energy (eV)','fontsize',15)
caxis([-5 0]+13)
caxis([-5 0]+13.5)
caxis([-5 0]+12-1)
% cblh2 = colorbar_labeled('Log_{10}(Flux(m^{-2}eV^{-1}s^{-1}))','log','fontsize',14);
cblh2 = colorbar_labeled('Flux (m^{-2}eV^{-1}s^{-1})','log','fontsize',14);
set(cblh2,'position',get(cblh2,'position')+[-0.01 0 0 0])
set(cblh2,'position',get(cblh2,'position')+[-0.01 0 0 0])
set(cblh2,'position',get(cblh2,'position')+[-0.01 0 0 0])
set(cblh2,'position',get(cblh2,'position')+[-0.01 0 0 0])
set(cblh2,'position',get(cblh2,'position')+[-0.01 0 0 0])
set(gca,'position',get(gca,'position')+[-0.01 0 0 0])
set(gca,'position',get(gca,'position')+[-0.01 0 0 0])
set(cblh2,'ytick',[1e6 1e8 1e10 1e12])
sp212pos = get(gca,'position');
set(gca,'position',[sp212pos(1:3),sp211pos(4)])
cb1h1pos = get(cblh2,'position')
set(cblh2,'position',[cb1h1pos(1:3),sp211pos(4)])
set(gca,'ytick',[10 100 1000])

%print -depsc2 -painters PnL-fig5n6-01.eps

axes(sp211)
axis([0 0.7 10 3002])
axes(sp212)
axis([0 0.7 10 3002])
% print -depsc2 -painters ../Figures/PnL-fig5n6-02.eps

%   Copyright © 2019 Bjorn Gustavsson, <bjorn.gustavsson@irf.se>
%   This is free software, licensed under GNU GPL version 2 or later
