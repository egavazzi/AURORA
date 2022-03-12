function [lon,lat,alt,xB,yB,zB] = getB_latlong2eq(lat0,long0,time_vec)
% getB_latlong2eq - magnetic field lines from Ramfjord using the IGRF.
% 
% getB_latlong2eq get the latitude, longitude altitude and the
% ECEF-cartesian coordinates of a magnetic-field-line from
% Ramfjordmoen, Tromsø, or any other location on Earth. 
% 
% Calling:
%  [lon,lat,alt,xB,yB,zB] = getB_latlong2eq(lat0,long0,time_vec)
% Input:
%  lat0     - latitude (deg), double scalar, positive northern
%             hemisphere, optional input, if empty defaults to
%             67.568 (location of Ramfjordmoen)
%  long0    - longitude (deg), double scalar, positive East of
%             Greenwitch, optional input, if empty defaults to
%             19.227 (location of Ramfjordmoen)
%  time_vec - optional input, defaults to [2007 7 17 6 30]
% Output:
%  lon - longitude (deg), double array 
%  lat - longitude (deg), double array 
%  alt - altitude above Earth surface (km), double array 
%  xB  - ECEF x-coordinate (km)
%  yB  - ECEF y-coordinate (km)
%  zB  - ECEF z-coordinate (km)
% All outputs have the same dimensions

%   Copyright © 2020 Bjorn Gustavsson, <bjorn.gustavsson@uit.no>
%   This is free software, licensed under GNU GPL version 2 or later


if nargin < 3 || isempty(time_vec)
  time = datenum([2007 7 17 6 30 0]);
else
  time = datenum(time_vec);
end

if isempty(lat0)
  lat_start = 67.586; % Geodetic latitudes in degrees.
else
  lat_start = lat0;
end
if isempty(long0)
  lon_start = 19.227; % Geodetic longitudes in degrees.
else
  lon_start = long0; % Geodetic longitudes in degrees.
end

alt_start = 0; % Altitude in km.

distance =  -83e3;% -50e3; % km.
nsteps = abs(distance)/1;

% Get the magnetic field line points.
[lat, lon, alt] = igrfline(time, lat_start, lon_start, alt_start, ...
                           'geod', distance, nsteps);

lat = lat(alt > -1);
lon = lon(alt > -1);
alt = alt(alt > -1);
% disp([alt(1:3)',min(alt),max(alt)])
lon(lon > 180) = lon(lon > 180) - 360;

[~,ieq] = max(alt);
% [ieq,size(alt)]
lat = lat(1:ieq);
lon = lon(1:ieq);
alt = alt(1:ieq);
% disp([alt(1:3)',min(alt),max(alt)])
% [ieq,size(alt),size(lon)]

if nargout > 3
  [x, y, z] = geod2ecef(lat, lon, alt*1e3);
  x = x/1e3;
  y = y/1e3;
  z = z/1e3;
  xB  = x;
  yB  = y;
  zB  = z;
end
