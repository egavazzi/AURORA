function h_out = cline(x,y,z,c,cmap,crange)
% CLINE - plots a 3D curve (x,y,z) encoded with scalar color data (c)
% using the specified colormap (default=jet);
%
% Calling:
%  h = cline(x,y,z,c,colormap,crange);
% 
% DBE 09/03/02
% Modified 20140814/BG: Added a CRANGE variable to make it possible
% to set the range for the colour span to work comparable to caxis
% for other pseudocolor plot functions.

if nargin == 0  % Generate sample data...
  x = linspace(-10,10,101);
  y = 2*x.^2+3;
  z = sin(0.1*pi*x);
  c = exp(z);
  % w = z-min(z)+1;
  cmap = 'jet';
elseif nargin<4
  fprintf('Insufficient input arguments\n');
  return;
end
if nargin < 5 || isempty(cmap)
  cmap = 'jet';
end
if nargin < 6 || isempty(crange)
  crange = [min(c) max(c)];
end

if ischar(cmap)
  cmap = colormap(cmap);                      % Set colormap
end

yy = linspace(crange(1),crange(2),size(cmap,1));  % Generate range of color indices that map to cmap
cm  =  interp1(yy,cmap,min(crange(2),max(crange(1),c)),'pchip')';         % Find interpolated colorvalue
cm(cm>1) = 1;                               % Sometimes iterpolation gives values that are out of [0,1] range...
cm(cm<0) = 0;


% Lot line segment with appropriate color for each data pair...
for i1 = (length(z)-1):-1:1,
  h(i1) = line([x(i1) x(i1+1)],[y(i1) y(i1+1)],[z(i1) z(i1+1)],'color',cm(:,i1),'LineWidth',3);
end

caxis(crange)
if nargout
  h_out = h;
end

return
