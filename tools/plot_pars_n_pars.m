function ph = plot_pars_n_pars(pars_of_ItMm,Yvals,idx2plot,idx4x,idx4Yvals,stdY,axlims,tstr,xstr,ystr,ftsz,dispP,lstl,clrs)
%                                         1     2        3     4         5    6      7    8    9   10   11,   12,
% PLOT_PARS_N_PARS - 
%   

%% ARGIN-checking section
if nargin < 5
  error('Too few input arguments')
end
if nargin < 6 || isempty(stdY)
  plot_w_plot = 1;
  % plot_w_errorbar = 0;
else
  % plot_w_errorbar = 1;
  plot_w_plot = 0;
end
if nargin < 7 || isempty(axlims)
  axlims = [];
end
if nargin < 8 || isempty(tstr)
  tstr = '';
end
if nargin < 9 || isempty(xstr)
  xstr = '';
end
if nargin < 10 || isempty(ystr)
  ystr = '';
end
if nargin < 11 || isempty(ftsz)
  ftsz = 14;
end
if nargin < 12 || isempty(dispP)
  dispP = 0;
end
if nargin < 13 || isempty(lstl)
  lstl = '.-';
end
if nargin < 13 || isempty(clrs)
  clrs{1} = 'r';
end

% keyboard

[~,idx4p] = sort(pars_of_ItMm(idx2plot,idx4x));
if plot_w_plot
  PH = plot(pars_of_ItMm(idx2plot(idx4p),idx4x),...
            Yvals(idx2plot(idx4p),idx4Yvals),...
            lstl,'linewidth',2,'markersize',18);
else
  PH = errorbar(repmat(pars_of_ItMm(idx2plot(idx4p),idx4x),size(idx4Yvals)),...
                Yvals(idx2plot(idx4p),idx4Yvals),...
                stdY(idx2plot(idx4p),idx4Yvals),...
                lstl,'linewidth',2,'markersize',18);
  %disp('Error-bars')
end
if ~isempty(tstr)
  title(tstr,'fontsize',ftsz)
end
if ~isempty(axlims)
  axis(axlims)
end
if ~isempty(ystr)
  ylabel(ystr,'fontsize',ftsz)
end
if ~isempty(xstr)
  xlabel(xstr,'fontsize',ftsz)
end
for iP = 1:numel(PH)
  set(PH(iP),'color',clrs{min(iP,end)})
end
grid on
if dispP ~= 0
  disp(tstr)
  disp(xstr)
  disp(pars_of_ItMm(idx2plot(idx4p),:))
  disp(ystr)
  disp(Yvals(idx2plot(idx4p),:))
  if dispP == -1
    disp('push any button')
    pause
  else
    pause(dispP)
  end
end

if nargout
  ph = PH;
end
