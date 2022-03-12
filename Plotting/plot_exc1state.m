function [ ph ] = plot_exc1state(t,h,q,title_str)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
subplot(2,1,1)
pcolor(t,h/1e3,q),shading flat
axis([t(1) t(end) 90 200])
title(title_str,'fontsize',15)
ylabel('height (km)','fontsize',15)
subplot(2,1,2)
ph = plot(t,sum(q),'linewidth',2);
xlabel('time','fontsize',15)
