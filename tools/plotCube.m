function H = plotCube(r_corners,varargin)
% PLOTCUBE - plot a cube/cuboid/rectangular parallelepiped
%   
% Calling:
%   H = plotCube(r_corners)


%   Copyright © 2019 Bjorn Gustavsson, <bjorn.gustavsson@irf.se>
%   This is free software, licensed under GNU GPL version 2 or later

h1 = plot3(r_corners([1:4,1,5:end,5],1),r_corners([1:4,1,5:end,5],2),r_corners([1:4,1,5:end,5],3),varargin{:});
h2 = plot3(r_corners([2,6],1),r_corners([2,6],2),r_corners([2,6],3),varargin{:});
h3 = plot3(r_corners([2,6]+1,1),r_corners([2,6]+1,2),r_corners([2,6]+1,3),varargin{:});
h4 = plot3(r_corners([2,6]+2,1),r_corners([2,6]+2,2),r_corners([2,6]+2,3),varargin{:});

if nargout
  H = [h1;h2;h3;h4];
end

