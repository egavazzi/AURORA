function r_a = r_arc(e1,e2,r0,r)
% R_ARC - points along a circular arc in 3-D
%   
% Calling:
%  r_a = r_arc(e1,e2,r0,r)
% Input:
%  e1 - unit-vector first end of arc, double array [1 x 3] or [3 x 1]
%  e2 - unit-vector second end of arc, double array [1 x 3] or [3 x 1]
%  r0 - centre-point or arc, double array [1 x 3] or  or [3 x 1],
%       optional argument, defaults to [0 0 0] 
%  r  - radius of arc, double scalar, optional argument, defaults
%       to one.
% Output:
%  r_a - points along arc, double array [n_points x 3], points are
%        separated by approximately 1 degree.
%
% Example:
%  e1 = [cos(12*pi/180)*cos(76*pi/180),sin(12*pi/180)*cos(76*pi/180),sin(76*pi/180)]
%  e2 = [cos(137*pi/180)*cos(-37*pi/180),sin(137*pi/180)*cos(-37*pi/180),sin(-37*pi/180)]
%  r0 = [1,1,1];
%  r  = 0.76;
%  R_arc = r_arc(e1',e2',r0,r);
%  plot3(r0(1),r0(2),r0(3),'r.','markersize',23)
%  hold on
%  plot3(R_arc(:,1),R_arc(:,2),R_arc(:,3),'k-','linewidth',2)


%  Copyright © Bjorn Gustavsson 20200304, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later


if nargin < 3 || isempty(r0)
  r0 = [0 0 0];
end
if nargin < 4 || isempty(r)
  r = 1;
end

e1 = e1(:)'/norm(e1);
e2 = e2(:)'/norm(e2);
e_to_rot_around = cross(e1,e2);

Theta = atan2(norm(cross(e1, e2)), dot(e1, e2));

n_points = round(Theta*180/pi)+1;
theta = linspace(0,Theta,n_points);


for i1 = n_points:-1:1,
  r_a(i1,:) = r0(:) + r*rot_around_v(e_to_rot_around,theta(i1))*e1';
end

