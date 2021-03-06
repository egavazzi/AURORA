clf
cx_lims = [-6 0] + 13;
subplot(1,1,1)                         
caxis(cx_lims)
colorbar_labeled('e^{-1}/m^2/s/eV/ster','log','interpreter','tex')
for i1 = 1:10,
  subplot(sph(i1,1),sph(i1,2),sph(i1,3)) 
  pcolor(t,E(1:size(Ie_top,3)),log10(squeeze(Ie_top(i1,:,:))'))
  shading flat
  cx(i1,:) = caxis;
  caxis(cx_lims)
  if i1 == 1
    title(['\theta_B: ',theta_strs{i1}],'interpreter','tex')
  else                   
    title([theta_strs{i1}],'interpreter','tex')             
  end 
  set(gca,...
      'YScale','log',...
      'TickDir','out',...
      'TickLength',[0.02 0.08],...
      'ytick',[10 100 1000])
end



for i1 = 6:10,
  subplot(sph(i1,1),sph(i1,2),sph(i1,3))
  xlabel('time (s)')
end
for i1 = [1 10],
  subplot(sph(i1,1),sph(i1,2),sph(i1,3)) 
  ylabel('energy (eV)')
end
