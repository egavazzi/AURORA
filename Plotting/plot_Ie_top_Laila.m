ftsz_ax = 12;
ftsz_xy = 14;
ftsz_t = 15;
dx_cbl = -0.01;
sx_cbl = -0.01;

idx4E_pa = [200 130 90];
if 1
  subplot(4,1,1)
  pcolor(t,E(1:size(Ie_top,3)),log10(squeeze(Ie_top(1,:,:))')),shading flat,set(gca,'yscale','log'),
  caxis([-6 0]+max(caxis))
  set(gca,'ytick',[10,30,100,300,1000])
  title('electron-flux || B','fontsize',ftsz_t)
  ylabel('energy (eV)','fontsize',ftsz_xy)
  set(gca,'xticklabel','','tickdir','out','fontsize',ftsz_ax,'box','off')
  axis([0 1.6 50 1675])
  
  subplot(4,1,2)
  pcolor(t,acos(mu_lims)*180/pi,log10(squeeze(Ie_top([1:end,end],:,idx4E_pa(1))))),shading flat
  title(sprintf('electron-flux vs pitch-angle E: %3.0f',E(idx4E_pa(1))),'fontsize',ftsz_t)
  caxis([-6 0]+max(caxis))
  ylabel('pitch-angle (deg)','fontsize',ftsz_xy)
  set(gca,'xticklabel','','tickdir','out','fontsize',ftsz_ax,'box','off')
  
  subplot(4,1,3)
  pcolor(t,acos(mu_lims)*180/pi,log10(squeeze(Ie_top([1:end,end],:,idx4E_pa(2))))),shading flat
  caxis([-6 0]+max(caxis))
  set(gca,'xticklabel','','tickdir','out','fontsize',ftsz_ax,'box','off')
  ylabel('pitch-angle (deg)','fontsize',ftsz_xy)
  title(sprintf('electron-flux vs pitch-angle E: %3.0f eV',E(idx4E_pa(2))),'fontsize',ftsz_t)

  subplot(4,1,4)
  pcolor(t,acos(mu_lims)*180/pi,log10(squeeze(Ie_top([1:end,end],:,idx4E_pa(3))))),shading flat,caxis([-6 0]+max(caxis))
  caxis([-6 0]+max(caxis))
  ylabel('pitch-angle (deg)','fontsize',ftsz_xy)
  title(sprintf('electron-flux vs pitch-angle E: %3.0f eV',E(idx4E_pa(3))),'fontsize',ftsz_t)
  xlabel('time (s)','fontsize',ftsz_xy)
  set(gca,'tickdir','out','fontsize',ftsz_ax,'box','off')
  drawnow
  subplot(4,1,1)
  cblh1 = colorbar_labeled('e^-1/eV/ster/m^2/s','log');
  set(cblh1,'position',get(cblh1,'position')+[dx_cbl 0 sx_cbl 0])
  set(cblh1,'ytick',[ 1e5 1e7 1e9])
  subplot(4,1,2)
  cblh2 = colorbar_labeled('e^-1/eV/ster/m^2/s','log');
  set(cblh2,'position',get(cblh2,'position')+[dx_cbl 0 sx_cbl 0])
  subplot(4,1,3)
  cblh3 = colorbar_labeled('e^-1/eV/ster/m^2/s','log');
  set(cblh3,'position',get(cblh3,'position')+[dx_cbl 0 sx_cbl 0])
  subplot(4,1,4)
  cblh4 = colorbar_labeled('e^-1/eV/ster/m^2/s','log');
  set(cblh4,'position',get(cblh4,'position')+[dx_cbl 0 sx_cbl 0])
end


set(cblh1,'ytick',[ 1e5 1e7 1e9])
set(cblh2,'ytick',[1e3 1e5 1e7]) 
set(cblh3,'ytick',[1e3 1e5 1e7])
set(cblh4,'ytick',[1e5 1e7 1e9])
