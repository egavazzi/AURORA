%% Phase-function-plots
% Script requires an energy-grid, suitably produced with the
% set-up-scripts

if ~exist('phfcnE_N2','var')
  [phfcnE_N2,phfcnI_N2] = phase_fcn_N2(E,(0:180)'*pi/180);
  [phfcnE_O2,phfcnI_O2] = phase_fcn_O2(E,(0:180)'*pi/180);
  [phfcnE_O,phfcnI_O] = phase_fcn_O(E,(0:180)'*pi/180);
end

Tftsz = 15;
Lftsz = 14;

sth = suptitle('Phase-functions');
set(sth,'fontsize',Tftsz)

%% N2
sph = subplot(3,2,1);
pcolor(E,(0:180),log10(phfcnE_N2)),shading flat
set(gca,'ytick',[0:30:180])
set(gca,'xscale','log','tickdir','out')
axis([5 1e4 0 180])
title('N_2 elastic collisions','fontsize',Tftsz)
ylabel('scattering angle (^\circ)','fontsize',Lftsz)
cblh = colorbar_labeled('','log');

sph(2) = subplot(3,2,2);
pcolor(E,(0:180),log10(phfcnI_N2)),shading flat
set(gca,'ytick',[0:30:180])
set(gca,'xscale','log','tickdir','out')
title('N_2 inelastic collisions','fontsize',Tftsz)
cblh(2) = colorbar_labeled('','log');

%% O2
sph(3) = subplot(3,2,3);
pcolor(E,(0:180),log10(phfcnE_O2)),shading flat
set(gca,'ytick',[0:30:180])
set(gca,'xscale','log','tickdir','out')
axis([5 1e4 0 180])
title('O_2 elastic collisions','fontsize',Tftsz)
ylabel('scattering angle (^\circ)','fontsize',Lftsz)
cblh(3) = colorbar_labeled('','log');

sph(4) = subplot(3,2,4);
pcolor(E,(0:180),log10(phfcnI_O2)),shading flat
set(gca,'ytick',[0:30:180])
set(gca,'xscale','log','tickdir','out')
title('O_2 inelastic collisions','fontsize',Tftsz)
cblh(4) = colorbar_labeled('','log');

%% O
sph(5) = subplot(3,2,5);
pcolor(E,(0:180),log10(phfcnE_O)),shading flat
set(gca,'ytick',[0:30:180])
set(gca,'xscale','log','tickdir','out')
axis([4 1e4 0 180])
title('O elastic collisions','fontsize',Tftsz)
ylabel('scattering angle (^\circ)','fontsize',Lftsz)
xlabel('Energy (eV)','fontsize',Lftsz)
cblh(5) = colorbar_labeled('','log');

sph(6) = subplot(3,2,6);
pcolor(E,(0:180),log10(phfcnI_O)),shading flat
set(gca,'ytick',[0:30:180])
set(gca,'xscale','log','tickdir','out')
title('O inelastic collisions','fontsize',Tftsz)
xlabel('Energy (eV)','fontsize',Lftsz)
cblh(6) = colorbar_labeled('','log');

drawnow

for i1 = 1:6,
  axpos = get(sph(i1),'position');
  cblp = get(cblh(i1),'position');
  set(cblh(i1),'position',[cblp(1:3),axpos(4)])
end
