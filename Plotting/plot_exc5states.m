function [ ph ] = plot_exc5states(t,Qall,title_str)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
for i1 = 5:-1:1,
  subplot(3,2,i1)
  q{i1} = sum(Qall{i1});
  ph = plot(t,q{i1},'linewidth',2);
  if i1 == 5
    xlabel('time','fontsize',15)
  end
  title(title_str{i1},'fontsize',15)
end
subplot(3,2,6)
for i1 = 5:-1:1,
  ph(i1) = plot(t,q{i1}/max(q{i1}),'color',rand(1,3),'linewidth',2);
  hold on
end
title('relative intensity')
xlabel('time','fontsize',15)
