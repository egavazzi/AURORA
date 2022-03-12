function [dtavg,stddt,dT] = dt_avgstd(t_ext1,t_ext2)
% DT_AVGSTD - average (and std) of difference between times in unsynchronised
%   sequences, of for example times for extremal values of time-varying
%   intensity values. A use-case is given two varying time-series
%   with local maxima and minima from a limited period where the
%   two time-series does not have the same number of local maxima
%   (minima), one time-series might be shifted just so that the
%   first, or last, maxima appears just outside the time-period of
%   observation. Then this function lines up the two
%   non-equal-length sequences with extreme-value-times giving the
%   smallest absolute average shift.
% 
% Calling:
%   [dtavg,stddt,dT] = dt_avgstd(t_ext1,t_ext2)
% Input:
%  t_ext1 - coordinates of extrema in sequence 1, double array
%           [1 x n1] or [n1 x 1] 
%  t_ext2 - coordinates of extrema in sequence 1, double array
%           [1 x n2] or [n2 x 1]
% Output:
%  dtavg - best average difference between the extrema in the sequences
%  stddt - standard deviation of the best difference between the
%          extrema in the sequences 
%  dT    - array with the differences between the best shifted
%          sequences 
%  t_ext1 and t_ext2 both have to be either column or row-vectors,
%  otherwise the new matlab "outer-difference" and "outer-sum"
%  behaviour might lead to unintended behaviour.
% 
% Example:
%   figure
%   w = 5;
%   t = linspace(0,5,501);           % time-array
%   y1 = erf(t) + 0.1*cos(w*t-0.5);  % time-series #1
%   y2 = erf(t) + 0.11*cos(w*t+0.5); % time-series #2
%   ph = plot(t,[y1;y2]);
%   set(ph(1),'color','b'),set(ph(2),'color','r')
%   k1 = local_max(y1);              % local maxima time-series #1
%   k1(k1==1|k1==numel(y1)) = [];    % but we cut local maxima at
%                                    % first and last point
%   k2 = local_max(y2);              % local maxima time-series #2
%   k2(k2==1|k2==numel(y2)) = [];    % and we cut local maxima at
%                                    % first and last point here too
%   hold on
%   plot(t(k1),y1(k1),'b.',t(k2),y2(k2),'r.','markersize',15)
%   [dtavg,stddt,dT] = dt_avgstd(t(k1),t(k2));
%   plot(t-dtavg,y2,'r--')
%   [dtavg2,stddt2,dT2] = dt_avgstd(t(k1),t(k2)-dtavg);
%   disp([dT;dT2]) 
%
% No argument checks or error-controls, if you want that you have
% to pay me good money.

% Copyright © 20190506 B. Gustavsson, <bjorn.gustavsson@uit.no>
%   This is free software, GPL version 3 or later applies

dtavg = nan;
stddt = nan;
dT{1} = [];

if isempty(t_ext1) || isempty(t_ext2)
  dtavg = nan;
  stddt = nan;
  dT    = [];
else
  n1 = numel(t_ext1);
  n2 = numel(t_ext2);
  dn = n2 - n1;
  if dn == 0
    dT    =     -(t_ext1 - t_ext2);
    dtavg = -mean(t_ext1 - t_ext2);
    stddt =   std(t_ext1 - t_ext2);
  elseif dn < 0
    dtLevel = inf;
    for i_shift = 0:abs(dn)
      dt_test = abs(mean(t_ext1((1:n2)+i_shift) - t_ext2));
      if dt_test < dtLevel
        dT      =     -(t_ext1((1:n2)+i_shift) - t_ext2);
        dtavg   = -mean(t_ext1((1:n2)+i_shift) - t_ext2);
        stddt   =   std(t_ext1((1:n2)+i_shift) - t_ext2);
        dtLevel = dt_test;
      end
    end
  else
    dtLevel = inf;
    for i_shift = 0:abs(dn)
      dt_test = abs(mean(t_ext1 - t_ext2((1:n1)+i_shift)));
      if dt_test < dtLevel
        dT      =     -(t_ext1 - t_ext2((1:n1)+i_shift));
        dtavg   = -mean(t_ext1 - t_ext2((1:n1)+i_shift));
        stddt   =   std(t_ext1 - t_ext2((1:n1)+i_shift));
        dtLevel = dt_test;
      end
    end
  end
end
