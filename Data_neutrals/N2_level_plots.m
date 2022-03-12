%% Script to plot a few of the N2 and N2+ energy levels

clf
ftsz = 14;
fntw = 'normal';% fntw = 'demi';
%% Plot the N_2 energy-levels
plot([0 3],[0 0],'k-','linewidth',2)
hold on
plot([5 8],[0 0],'k-','linewidth',2)
plot([6 8],N2_levels(13,[1 1]),'k-','linewidth',2)
plot([6 8],N2_levels(14,[1 1]),'k-','linewidth',2)
plot([6 8],N2_levels(22,[1 1]),'k-','linewidth',2)
plot([0 2],N2_levels(21,[1 1]),'k-','linewidth',2)
plot([0 2],[1 1]*12.5,'k-','linewidth',2)
%% Plot the N_2^+ energy-levels
plot([3 5],N2_levels(end-2,[1 1]),'k-','linewidth',2)  
plot([3 5],N2_levels(end-1,[1 1]),'k-','linewidth',2)
plot([3 5],N2_levels(end,[1 1]),'k-','linewidth',2)
%% Label the energy-levels
th = text(3.75,0,'N_2(X^1\Sigma_g^+)','fontsize',ftsz,'fontweight',fntw);
th(end+1) = text(4.75,N2_levels(13,1)-.25,'N_2(A^3\Sigma_u^+)',...
                 'fontsize',ftsz,'fontweight',fntw);
th(end+1) = text(4.75,N2_levels(14,1)-.25,'N_2(B^3\Pi_g)',...
                 'fontsize',ftsz,'fontweight',fntw);
th(end+1) = text(4.75,N2_levels(22,1)-.25,'N_2(C^3\Pi_u)',...
                 'fontsize',ftsz,'fontweight',fntw);
th(end+1) = text(2.25,N2_levels(21,1)-.25,'N_2(a^1\Pi_g)',...
                 'fontsize',ftsz,'fontweight',fntw);
th(end+1) = text(2.25,12.5-.25,'N_2(b^1\Pi_u)',...
                 'fontsize',ftsz,'fontweight',fntw);
th(end+1) = text(5.25,N2_levels(end-2,1)-.25,'N_2^+(X^2\Sigma_g^+)',...
                 'fontsize',ftsz,'fontweight',fntw);
th(end+1) = text(5.25,N2_levels(end-1,1)-.25,'N_2^+(A^2\Pi_u)',...
                 'fontsize',ftsz,'fontweight',fntw);
th(end+1) = text(5.25,N2_levels(end,1)-.25,'N_2^+(B^2\Sigma_u^+)',...
                 'fontsize',ftsz,'fontweight',fntw);
%% Arrow the prominent band-transitions
% For matlab-versions where arrow works:
try
  arrow([7 N2_levels(13,1)],[7,0])
  arrow([6.5 N2_levels(14,1)],[6.5 N2_levels(13,1)],'facecolor','r','edgecolor','r')
  arrow([7.5 N2_levels(22,1)],[7.5 N2_levels(14,1)],'facecolor',[.5 0 .8],'edgecolor',[.5 0 .8])
  arrow([.75 N2_levels(21,1)],[.75 0])
  arrow([.25,12.5],[.25 N2_levels(21,1)])
  arrow([1.5,12.5],[1.5 0]) 
  arrow([3.25 N2_levels(end-1,1)],[3.25 N2_levels(end-2,1)],'facecolor','r','edgecolor','r')
  arrow([4.75 N2_levels(end,1)],[4.75 N2_levels(end-2,1)],'facecolor',[.2 0 1],'edgecolor',[.2 0 1])
catch
  % for versions where arrow doesn't work try with arrow3
  ah1 = arrow3([7 N2_levels(13,1)],[7,0],'k-',.8,1.25);
  ah2 = arrow3([6.5 N2_levels(14,1)],[6.5 N2_levels(13,1)],'r-',.8,1.25);
  ah3 = arrow3([7.5 N2_levels(22,1)],[7.5 N2_levels(14,1)],'^i-',.8,1.25);
  ah4 = arrow3([.75 N2_levels(21,1)],[.75 0],'k-',.8,1.25);
  ah5 = arrow3([.25,12.5],[.25 N2_levels(21,1)],'k-',.8,1.25);
  ah6 = arrow3([1.5,12.5],[1.5 0],'k-',.8,1.3);
  ah7 = arrow3([3.25 N2_levels(end-1,1)],[3.25 N2_levels(end-2,1)],'r-',.8,1.25);
  ah8 = arrow3([4.75 N2_levels(end,1)],[4.75 N2_levels(end-2,1)],'^i-',.8,1.25);
  set(ah1(1),'color','k'),set(ah1(2),'facecolor','k')
  set(ah2(1),'color','r'),set(ah2(2),'facecolor','r')
  set(ah3(1),'color',[.6 .1 1]),set(ah3(2),'facecolor',[.6 .1 1])
  set(ah4(1),'color','k'),set(ah4(2),'facecolor','k')
  set(ah5(1),'color','k'),set(ah5(2),'facecolor','k')
  set(ah6(1),'color','k'),set(ah6(2),'facecolor','k')
  set(ah7(1),'color','r'),set(ah7(2),'facecolor','r')
  set(ah8(1),'color',[.4 0 1]),set(ah8(2),'facecolor',[.4 0 1])
end
  %% Label the bands
text(6.6183,1.5018,'Vegard-Kaplan','rotation',90,...
     'fontsize',ftsz,'fontweight',fntw)
text(6.65,6.6835,'1st positive',...
     'fontsize',ftsz,'fontweight',fntw)
text(6.0651,9.1414,'2nd positive',...
     'fontsize',ftsz,'fontweight',fntw)
text(2.0911,16.0101,'Meinel',...
     'fontsize',ftsz,'fontweight',fntw)
text(3.25,17.569,'1st negative',...
     'fontsize',ftsz,'fontweight',fntw)
text(0.4710,1.9411,'Lyman-Birge-Hopfield','rotation',90,...
     'fontsize',ftsz,'fontweight',fntw)
text(1.7920,3.7414,'Birge-Hopfield','rotation',90,...
     'fontsize',ftsz,'fontweight',fntw)
text(0.5640,9.8419,'Janin','rotation',90,...
     'fontsize',ftsz,'fontweight',fntw)
%% Title-n-label the diagram
title('N_2 and N_2^+ energy-levels','fontsize',15)
ylabel('Energy (eV)','fontsize',14)
set(gca,'xtick',[],'box','on','ytick',0:20)
axis([0 8 -.5 20])
