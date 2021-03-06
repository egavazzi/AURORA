figure('position',[2, 554, 1238, 420])
c_o_mu = mu_avg(mu_lims);
spp = [2*ones(numel(BeamW),1),...
       numel(BeamW)/2*ones(numel(BeamW),1),...
       [1:numel(BeamW)/2,numel(BeamW):-1:(numel(BeamW)/2+1)]'];
for i1 = 1:numel(BeamW)
  subplot(spp(i1,1),spp(i1,2),spp(i1,3))
  imagesc(0:180,180:-1:0,Pmu2mup(:,:,i1))
  title(sprintf('%3.1f - %3.1f',180-theta_lims(i1),180-theta_lims(i1+1)))
end
colormap(inferno)
figure
subplot(1,2,1)
imagesc(squeeze(sum(Pmu2mup,3))),colorbar
subplot(1,2,2)
imagesc(squeeze(sum(Pmu2mup,3))-1),colorbar
figure
subplot(3,1,1)
imagesc(0:180,1:18,2*pi*theta2beamW)
title('Beam-weights')
subplot(3,1,2)
plot(BeamW*2*pi,'s-')
title('Beam-widths (ster)')
grid on
subplot(3,1,3)
plot(c_o_mu,'s-')
title('\mu-avg')
xlabel('stream number')
grid on