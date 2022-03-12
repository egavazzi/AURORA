function [ ph ] = plot_exc5states(t,Qall,title_str,spp)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
for i1 = numel(Qall):-1:1,
  subplot(spp(i1,1),spp(i1,2),spp(i1,3))
  q{i1} = sum(Qall{i1});
  ph = plot(t,q{i1},'linewidth',2);
  if i1 == numel(Qall)
    xlabel('time','fontsize',15)
  end
  title(title_str{i1},'fontsize',15)
end
subplot(spp(end,1),spp(end,2),spp(end,3))
for i1 = numel(Qall):-1:1,
  ph(i1) = plot(t,q{i1}/max(q{i1}),'color',rand(1,3),'linewidth',2);
  hold on
end
try
  cmlines(ph)
end
legend(ph,title_str)
title('relative intensity')
xlabel('time','fontsize',15)
