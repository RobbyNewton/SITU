% LoadEASE25
% Loads the binary file of lats and lons
% of the N.H. 25-km EASE grid. 
% The grid is set up for IDL, which 
% reads left to right for "x" and 
% bottom to top for "y".  I.e.: the 
% row and column indices are opposite of 
% the matlab (matrix) convention and the 
% rows are indexed in the opposite direction. 
% (Or ... IDL is "cartesian" whereas Matlab 
% is "matrix".)

addpath('/storage/patricia/ForecastingArcticIceTracks/EASEGrid')
fid = fopen('NLLATLSB','r') ;
NHLat25 = fread(fid,'int32') ;
fclose(fid) ;
NHLat25 = NHLat25/100000 ;
NHLat25(NHLat25==NHLat25(1)) = NaN ;
NHLat25 = reshape(NHLat25, 721, 721) ;

fid = fopen('NLLONLSB','r') ;
NHLon25 = fread(fid,'int32') ;
fclose(fid) ;
NHLon25 = NHLon25/100000 ;
NHLon25(NHLon25==NHLon25(1)) = NaN ;
NHLon25 = reshape(NHLon25, 721, 721)' ;
