%% Script generating t-z plots of differential electron-fluxes


try
  [t_PnL,h_atm,E,mu_lims,IeZTE,mu_scatterings] = Ie_ztE_loader({'PnLs-05-4000-3keV-0-10'});
catch
  [t_PnL,h_atm,E,mu_lims,IeZTE,mu_scatterings] = Ie_ztE_loader({'PnLs-05-4000-3keV-FA'});
end

dE = diff(E);dE(end+1) = dE(end);
BeamW = 2*pi*mu_scatterings{3};

for i1 = numel(BeamW):-1:1,
  theta_str{i1} = sprintf('%3.1f - %3.1f',...
                          180-180/pi*acos(mu_lims(i1)),...
                          180-180/pi*acos(mu_lims(i1+1)));
end

theta_strs'

spp = [2*ones(numel(c_o_mu),1),5*ones(numel(c_o_mu),1),[(1:5)';(10:-1:6)']];

[dE90,iE90] = min(abs(E-90))    
[dE100,iE100] = min(abs(E-100))  
[dE500,iE500] = min(abs(E-500))  
[dE1000,iE1000] = min(abs(E-1000))
[dE2000,iE2000] = min(abs(E-2000))

E([iE90,iE100,iE500,iE1000,iE2000])

figure('position',[692, 556, 997, 418])
plot_IeztE_3DtzofE(t_PnL,h_atm/1e3,E,IeZTE,dE,BeamW*2*pi,[6.5 12],spp, theta_strs,iE100);
figure('position',[692, 556, 997, 418])
plot_IeztE_3DtzofE(t_PnL,h_atm/1e3,E,IeZTE,dE,BeamW*2*pi,[6.5 12],spp, theta_strs,iE90); 
figure('position',[692, 556, 997, 418])
plot_IeztE_3DtzofE(t_PnL,h_atm/1e3,E,IeZTE,dE,BeamW*2*pi,[7 12.5],spp, theta_strs,iE500);
figure('position',[692, 556, 997, 418])
plot_IeztE_3DtzofE(t_PnL,h_atm/1e3,E,IeZTE,dE,BeamW*2*pi,[7 12.5],spp, theta_strs,iE1000);
figure('position',[692, 556, 997, 418])
plot_IeztE_3DtzofE(t_PnL,h_atm/1e3,E,IeZTE,dE,BeamW*2*pi,[7 12.5],spp, theta_strs,iE2000);


%   Copyright © 2019 Bjorn Gustavsson, <bjorn.gustavsson@irf.se>
%   This is free software, licensed under GNU GPL version 2 or later
