function [B,lB,lon,lat,alt,xB,yB,zB] = make_BofL(lat0,long0,Btime,z_0,z_end)
% MAKE_BOFL - Magnetic field strength, magnetic field-line-length
% and geographical coordinates (GEO lat-long-alt, ECEF x-y-z)
% 
% Calling:
%  [B,lB,lon,lat,alt,xB,yB,zB] = make_BofL(lat0,long0,Btime,z_0,z_end)
% Input:
%  lat0,long0,Btime,z_0,z_end
%  lat0     - latitude (deg), double scalar, positive northern
%             hemisphere, optional input, if empty defaults to
%             67.568 (location of Ramfjordmoen)
%  long0    - longitude (deg), double scalar, positive East of
%             Greenwitch, optional input, if empty defaults to
%             19.227 (location of Ramfjordmoen)
%  Btime    - optional input, defaults to [2007 7 17 6 30]
%  z_0      - lowest altitude (m), double scalar
%  z_end    - hihgest altitude (m), double scalar
% Output:
%  B   - Magnetic-field-strength (T), double array
%  lB  - distance along magnetic field-line (m), double array.
%  lon - longitude (deg), double array 
%  lat - longitude (deg), double array 
%  alt - altitude above Earth surface (km), double array 
%  xB  - ECEF x-coordinate (km)
%  yB  - ECEF y-coordinate (km)
%  zB  - ECEF z-coordinate (km)
% All outputs have the same dimensions

%   Copyright © 2020 Bjorn Gustavsson, <bjorn.gustavsson@uit.no>
%   This is free software, licensed under GNU GPL version 2 or later

if nargin < 1 || isempty(lat0)
  latitude = 69.58;
else
  latitude = lat0;
end
if nargin < 2 || isempty(long0)
  longitude = 19.23;
else
  longitude = long0;
end
if nargin < 3 || isempty(Btime)
  t4B = [2005 08 10 0 0 0];
else
  t4B = Btime;
end

% $$$ [lon,lat,alt,xB,yB,zB] = plotB_latlong_conjugates(latitude,...
% $$$                                                   longitude,...
% $$$                                                   t4B);

[lon,lat,alt,xB,yB,zB] = getB_latlong2eq(latitude,longitude,t4B);
% disp(alt([1 2 3 end-2:end])')
[dz1,i_of_interest(1)] = min(abs(alt*1e3-z_0));
[dz2,i_of_interest(2)] = min(abs(alt*1e3-z_end));

lon = lon(i_of_interest(1):i_of_interest(2)); 
lat= lat(i_of_interest(1):i_of_interest(2));
alt= alt(i_of_interest(1):i_of_interest(2));
xB= xB(i_of_interest(1):i_of_interest(2));
yB= yB(i_of_interest(1):i_of_interest(2));
zB= zB(i_of_interest(1):i_of_interest(2));

[Bt, Bp, Br] = igrf(datenum(t4B(1:3)),lon,lat,alt);
B = (Bt.^2+Bp.^2+Br.^2).^(1/2)*1e-9;
% B = interp1(alt(1:800),B0(1:800),h_atm/1e3,'pchip');

rXYZ = [xB,yB,zB]*1e3;
dRB = (diff(rXYZ));

lB = sum(dRB.^2,2).^.5;
lB = [0;cumsum(lB)];
