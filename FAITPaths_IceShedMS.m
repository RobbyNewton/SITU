% FAITPaths_IceShedMS
%
% A script to set paths for the Sea Ice Trcking Utility scripts.
%
% This file will have to be edited to your specific environment.

% datapath = Where is the input data stored?

% sl = '\' ;  % for Windows. Set in startup.m
ResultPath = 'Where will the output be stored?' ;

EEZPath             = [dbpath 'LITS' sl 'EEZ' sl] ;
BathPath            = [dbpath 'LITS' sl 'Matfiles' sl 'EASEGrid' sl] ;
GridPath            = BathPath ;
LITSDatLocal        = [basepath 'LITSData' sl 'CESM-LE' sl] ; % Downloaded from polarapps.
LITSDatLocalLW      = [basepath 'LITSData' sl 'CESM-LW' sl] ; % Mirrored on polarapps.
CESMInterpPath      = [LITSDatLocal 'Interp' sl] ; % V3, the 2019 version.  
CESMOutputPath      = [LITSDatLocal 'EASEGrid' sl] ;
CESMLWInterpPath    = [LITSDatLocalLW 'Monthly' sl 'Interp' sl] ;
CESMLWOutputPath    = [LITSDatLocalLW 'Monthly' sl 'EASEGrid' sl] ;
CESMTrackedIce      = [LITSDatLocal 'TrackedIce' sl] ;
CESMLWTrackedIce    = [LITSDatLocalLW 'TrackedIce' sl] ;

addpath(BathPath, GridPath, LITSDatLocal, LITSDatLocalLW, matpaths.LITS) 
addpath([matpaths.LITS 'polygeom'])
addpath([matpaths.LITS 'PolygonClipper'])
addpath([dbpath 'LITSData' sl 'XYLagsV2' sl])        % matfiles with lagrangian tracks.  
% The PPF components, re-writen in .mat format.
PPFComponents   = [datpaths.LITS 'PolarPathfinder\RawDriftVectorsMATFormat'] ; 
