function [lats, lons] = XY2LatLon(xs,ys)
% XY2LatLon
% [lats, lons] = XY2LatLon(xs,ys)
% Translates from 25-km EASE grid coordinates xs, ys
% to latitude and longitude coordinates.  
%
% xs and ys can be scalars, vectors or matrices.
% 
% B Newton, Jan. 2018.

north_x_y_lat_lon = ...
load('C:\Users\home\Dropbox (LDEO)\LITS\Matfiles\EASEGrid\north_x_y_lat_lon') ;
% x       = north_x_y_lat_lon(:,1) + 1 ; % 0-based to 1-based.
% y       = north_x_y_lat_lon(:,2) + 1 ;
lat     = north_x_y_lat_lon(:,3) ;
lon     = north_x_y_lat_lon(:,4) ;
% Puth the coordinates into 2-D grid.
% XX  = reshape(x,361,361) ;
% YY  = reshape(y,361,361) ;
LAT = reshape(lat,361,361) ;
LON = reshape(lon,361,361) ;
clear north_x_y_lat_lon

lats    = interp2(LAT,ys,xs) ;
lons    = interp2(LON,ys,xs) ;