if 1
  subplot(4,1,1)
  pcolor(t,E(1:size(Ie_top,3)),log10(squeeze(Ie_top(1,:,:))')),shading flat,set(gca,'yscale','log'),
  caxis([-6 0]+max(caxis))
  set(gca,'ytick',[10,30,100,300,1000])
  title('electron-flux || B','fontsize',15)
  ylabel('energy (eV)','fontsize',15)
  cblh1 = colorbar_labeled('e^-1/eV/ster/m^2/s','log');
  set(cblh1,'position',get(cblh1,'position')+[-0.02 0 0 0])
  set(gca,'xticklabel','','tickdir','out')
  
  subplot(4,1,2)
  pcolor(t,acos(mu_lims)*180/pi,log10(squeeze(Ie_top([1:end,end],:,200)))),shading flat
  title(sprintf('electron-flux vs pitch-angle E: %3.1f',E(200)),'fontsize',15)
  ylabel('pitch-angle (deg)','fontsize',15)
  cblh2 = colorbar_labeled('e^-1/eV/ster/m^2/s','log');
  set(gca,'xticklabel','','tickdir','out')
  set(cblh2,'position',get(cblh2,'position')+[-0.02 0 0 0])
  
  subplot(4,1,3)
  pcolor(t,acos(mu_lims)*180/pi,log10(squeeze(Ie_top([1:end,end],:,130)))),shading flat
  caxis([-6 0]+max(caxis))
  set(gca,'xticklabel','','tickdir','out')
  ylabel('pitch-angle (deg)','fontsize',15)
  cblh3 = colorbar_labeled('e^-1/eV/ster/m^2/s','log');
  set(cblh3,'position',get(cblh3,'position')+[-0.02 0 0 0])
  title(sprintf('electron-flux vs pitch-angle E: %3.1f eV',E(130)),'fontsize',15)

  subplot(4,1,4)
  pcolor(t,acos(mu_lims)*180/pi,log10(squeeze(Ie_top([1:end,end],:,90)))),shading flat,caxis([-6 0]+max(caxis))
  caxis([-6 0]+max(caxis))
  ylabel('pitch-angle (deg)','fontsize',15)
  title(sprintf('electron-flux vs pitch-angle E: %3.1f eV',E(90)),'fontsize',15)
  xlabel('time (s)','fontsize',15)
  set(gca,'tickdir','out')
  cblh4 = colorbar_labeled('e^-1/eV/ster/m^2/s','log');
  set(cblh4,'position',get(cblh4,'position')+[-0.02 0 0 0])
end
