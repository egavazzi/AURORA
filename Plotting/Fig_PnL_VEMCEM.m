%% Reproduce Volume and column emission plot of Peticolas and Lummerzheim
% This script reproduces the Volume and column emissions plot in
% Figures 3 in Peticolas and Lummerzheim (2000), with fluxes
% from AURORA. This script requires the results from a 0.7 s long
% integration of 5 Hz square-wave modulated FAB from 4000 km. This
% is generated with one of the parameter-sets in
% Loop_finalRun_PnL_flickering10beams.m, that will be automatically
% saved in a directory named PnLs-05-4000-3keV-0-10


cd PnLs-05-4000-3keV-FA

load Qzt_all_L.mat
load I_lambda_of_t.mat
t_PnL = t;
cd ..
figure('position',[517 398 723 576])
subplot(2,1,2)
plot(t_PnL,I_4278/1e10,'linewidth',2)
set(gca,'fontsize',12)
xlabel('time (s)','fontsize',15)
ylabel('Intensity (R)','fontsize',15)
MAP4278 = (max(I_4278) - min(I_4278(end/2:end)))/max(I_4278);
tstr = sprintf('4278 A Column Emission Rate: %3.1f%% modulation',MAP4278*100);
title(tstr,'fontsize',16)
set(gca,'box','off','tickdir','out')
subplot(2,1,1)
pcolor(t_PnL,h_atm/1e3,Q4278),shading flat
set(gca,'fontsize',12)
axis([0 0.7 100 250])
set(gca,'box','off','tickdir','out','xticklabel','')
cblh1 = colorbar_labeled('Emission Rate(m^{-3}s^{-1})','linear','fontsize',14);

set(cblh1,'position',get(cblh1,'position')+[-0.01 0 0 0])
set(cblh1,'position',get(cblh1,'position')+[-0.01 0 0 0])
title('4278 A Volume Emission Rate','fontsize',16)
ylabel('Height (km)','fontsize',15)

delete(cblh1)
cblh1 = colorbar_labeled('Emission Rate(m^{-3}s^{-1})','linear','fontsize',14);
set(cblh1,'position',get(cblh1,'position')+[-0.01 0 0 0])
set(cblh1,'position',get(cblh1,'position')+[-0.01 0 0 0])
%  print -depsc2 -painters ../Figures/PnL-fig_VEMCEM-01.eps


%   Copyright © 2019 Bjorn Gustavsson, <bjorn.gustavsson@irf.se>
%   This is free software, licensed under GNU GPL version 2 or later
