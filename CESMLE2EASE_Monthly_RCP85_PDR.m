% CESMLE2EASE_Monthly_RCP85_PDR
%
% A script to:
%   * Load CESMLE RCP 8.5 model output
%   * Interpolate variables (SIC, SIT, SIT categories, SLP) onto EASE Grid
%   * Rotate the velocity vectors towards the real North Pole
%   * Interpolate onto EASE Grid
%   * Store the interpolants and fields in a set of monthly files.
%     This will pre-calculate for later use.
%
% In the ice concentration files:
% Time = time has 1140 months.
% Latitude = TLAT has 104 elements.
% Longitude = TLON has 320 elements.
% The data fields, e.g.:
% aice has 320 X 104 X 1140 elements = TLON x TLAT x time.
%
% In the velocity files:
% Time = time has 1140 months.
% Latitude = lat has 104 elements.
% Longitude = lon has 320 elements.
% The data fields, e.g.:
% uvel has 320 X 104 X 1140 elements = lon x lat x time.
%
% In the sea-level pressure files:
% Time = time has 1140 months.
% Latitude = lat has 192 elements.
% Longitude = lon has 288 elements.
% The data fields, e.g.:
% PSL has 288 X 192 X 1140 elements = lon x lat x time. 
%
% The data are stored in their native grid.
% U is positive toward the East in the NATIVE grid
% with the North Pole over central northern Greenland.
% V is positive towards THAT pole.  
% In at least some files there is a discontinuity along the
% N. Alaskan coast ... that feature is in the data.
% And there is a strong offshore velocity along the SE Greenland Current.
% And we have a question about the flows off Iceland.  
%
% ULAT and ULON, which we use:
% latU = double(ncread(filenameU, 'ULAT',[1 17],[inf inf])) ;
% are the REAL, GEOPHYSICAL latitude and longitude at the 
% CCSM grid point.  
% 
% To plot the CCSM Output in its native grid:
%   figure
%   makearcticframeNoTopo
%   pcolorm(latU, lonU, vice), figure(gcf)
%   colorbar, axis equal
%   quiverm(latU, lonU, uice, vice), figure(gcf)
% and then the same for vice. 
%
% EASE GRID COORDINATES:
% The grid starts with (0,0) on the top left corner,
% and proceeds to the right and downward. (Top Left means
% top and left when the map is viewed so that the NP is
% in the center and the Grn.Merid. points downward ... i.e.:
% Top left is in the Pacific.
% Each time one loads the grid, add 1 to x and y so they
% can serve as indices in Matlab.
% The grid is 361 X 361 and centered on the N. Pole.
% The Greenwich Meridian is pointing from the center down.
% x increases from left to right.  It is a natural counter of
% rows in IDL or columns in Matlab.
% y increases from top to bottom.
% u is positive from left to right.  (+x direction!, N. Amer to Eurasia.)
% v is positive from bottom to top.  (-y direction!, Atlantic to Pacific.)
%
% C-Grid. 
%
% For the sea-ice concentration files, TLAT and TLON have NaN values where
% there is a NaN in aice (over the continents). This was not the case when
% working with the CCSM data, as the lat and lon matrices are complete.
% This brings a problem when creating the interpolant function as NaN
% values of lat/lon are not permitted. Since the TLAT/TLON matrices are
% exactly the same for CESM and CCSM except over the continents, we will
% use the ones from CCSM to load the lat/lon (same for ULAT and ULON with
% the ice velocities).
%
% For the thickness categories, the data is stored in the same way as
% sea-ice concentration. It reprensents the percentage of ice concentration
% in each thickness categories, the sum of all categories adding up to the
% full sea-ice concentration field. More information can be found in the
% file cicedoc.pdf
%
% Patricia, July 2015
% Before saving the Interpolants for u- and v- component of the velocity,
% we set all NaN values to zeros. This is because of the difference in the 
% way the grid is saved between the vector fields and the scalar field. On
% a C-Grid, we would think that the grid knot of the continents (NaNs) in 
% the sic field (scalar) would correspond to a zero (boundary condition) in
% the u- and v- field (vector) and that the continents (NaNs) would be
% slightly smaller in the vector field. This is not the case. The scalar
% fields and the vector fields have the same grid size (we would think u-
% and v- should have one more row and column). In the u- and v- fields, the
% continents (NaNs) are the same as the ones in the sic field except that
% they are shifted 1/2 of a grid point towards the North Pole of the grid
% (over Greenland). The maximum latitude in the vector fields is slightly
% greater than the one in the scalar fields. On the downside of the
% continents (away from the North Pole of the grid), velocities are always
% zero. This created a discrepancy between the sic field and the u- and v-
% fields when tracking the ice and using the Interpolant (some EASEGrid
% locations would have sic = NaN and u&v /= NaN or vice versa). To avoid
% this problem, we set all NaNs in the u- and v- fields to zero and
% consider the continents (NaN values) in the sic field (scalar) as the
% true continents. For a more visual understanding, figures can be found in
% in the Dropbox under /ForecastingArcticIceTracks/Results_PDR/CCSM_Grid/
%
% Patricia, April 2015
% Last update: May 2020

%% Common variables
% % When working on Mac
% addpath('/Users/pade7652/Google Drive/CUBoulder/Projects/EEZFuture/Scripts')
% FAITPaths_Mac_PDR

% When working on Cheyenne
FAITPaths_Cheyenne

% Ensemble members
member = {'001' '002' '003' '004' '005' '006' '007' '008' '009' '010' ...
    '011' '012' '013' '014' '015' '016' '017' '018' '019' '020' ...
    '021' '022' '023' '024' '025' '026' '027' '028' '029' '030' ...
    '031' '032' '033' '034' '035' '101' '102' '103' '104' '105'} ;

% Read EASE Grid coordinates
load('north_x_y_lat_lon')
x       = north_x_y_lat_lon(:,1) + 1 ;
y       = north_x_y_lat_lon(:,2) + 1 ;
lat     = north_x_y_lat_lon(:,3) ;
lon     = north_x_y_lat_lon(:,4) ;

XX      = reshape(x,361,361) ;
YY      = reshape(y,361,361) ;
LAT     = reshape(lat,361,361) ;
LON     = reshape(lon,361,361) ;
l = 125 ;
clear north_x_y_lat_lon

% % Load CESM mask
% load('CESMmask.mat') % Ocean = 1; land = 0, per CESM grid
% 
% % New CESM mask for surface plot
% new_CESMmask = CESMmask ;
% new_CESMmask(new_CESMmask==1) = NaN ;
% new_CESMmask(new_CESMmask==0) = 1 ;
% 
% % Load bathymetry
% load('bathyEASE')

% % CCSM file for CCSM-Grid Lats and Lons
% filename_CCSM = ([CCSMRawPath 'RCP8.5/b40.rcp8_5.1deg.007.cice.h.aice_nh.200501-210012.nc']) ;
% 
% % CCSM-grid lats and lons !! -> see comment section above
% % for CCSM SIC files: all Lat = north of 30N (matches 361x361 grid)
% latT = double(ncread(filename_CCSM,'TLAT',[1 1],[inf inf])) ;
% lonT = double(ncread(filename_CCSM,'TLON',[1 1],[inf inf])) ;
% lonT(lonT>180) = lonT(lonT>180)-360 ;
% 
% latU = double(ncread(filename_CCSM,'ULAT',[1 1],[inf inf])) ;
% lonU = double(ncread(filename_CCSM,'ULON',[1 1],[inf inf])) ;
% lonU(lonU>180) = lonU(lonU>180)-360 ;
% 
% % Local angle from CCSM grid to geographical grid
% angleU = double(ncread(filename_CCSM,'ANGLE',[1 1],[inf inf])) ;

%% Looping through all the ensemble members
tic
for imem = 1:length(member)
    disp(['Member ' member{imem}])
    
%    % Read in the CESM data
%    filenameC = ([CESMLERawMPath 'RCP8.5/aice/b.e11.BRCP85C5CNBDRD.f09_g16.' ...
%        member{imem} '.cice.h.aice_nh.200601-210012.nc']) ;
%    filenameU = ([CESMLERawMPath 'RCP8.5/uvel/b.e11.BRCP85C5CNBDRD.f09_g16.' ...
%        member{imem} '.cice.h.uvel_nh.200601-210012.nc']) ;
%    filenameV = ([CESMLERawMPath 'RCP8.5/vvel/b.e11.BRCP85C5CNBDRD.f09_g16.' ...
%        member{imem} '.cice.h.vvel_nh.200601-210012.nc']) ;
%    
%    % Declare initial variables
%    sicEASE = NaN(361,361,12) ;
%    uiceEASE = NaN(361,361,12) ;
%    viceEASE = NaN(361,361,12) ;
%    uiceEASExy = NaN(361,361,12) ;
%    viceEASExy = NaN(361,361,12) ;
%    ystart = 2006 ; 
%    
%    %% Looping through the years
%    for iy = 2006:2100
%        disp(['Year: ' int2str(iy)])
%        
%        % Looping through the months
%        for im = 1:12
%            disp(['Year: ' int2str(iy) ', Month: ' int2str(im)])
%            mon = (iy-ystart)*12 + im ;
%            
%            %% Sea-ice concentration
%            % Read in the CESM data
%            sic = ncread(filenameC,'aice',[1 1 mon],[inf inf 1]) ;
%            
%            % Create interpolant function
%            Fc = scatteredInterpolant(latT(:), lonT(:), sic(:)) ;
%            
%            % Interpolate to the EASE Grid
%            sicEASE(:,:,im) = Fc(LAT,LON) ;
%            
%            % Discontinuity problem
%            sicEASE(181,1:180,im) = NaN ;
%            sicEASE(180,1:59,im) = NaN ;
%            sicEASE(182,1:59,im) = NaN ;
%            sicEASE(180,180,im) = NaN ;
%            sicEASE(182,181,im) = NaN ;
%            sicEASE(180,181,im) = NaN ;
%            sicEASE(181,181,im) = NaN ;
%            sicEASE(180,182,im) = NaN ;
%            sicEASE(181,182,im) = NaN ;
%            
%            sicEASE(180,1:59,im) = nanmean(sicEASE(178:183,1:59,im),1) ;
%            sicEASE(182,1:59,im) = nanmean(sicEASE(180:184,1:59,im),1) ;
%            sicEASE(181,1:180,im) = nanmean(sicEASE(180:182,1:180,im),1) ;
%            sicEASE(180,180,im) = nanmean(sicEASE(179:181,180,im),1) ;
%            sicEASE(182,181,im) = nanmean(sicEASE(182,180:182,im),2) ;
%            sicEASE(180,181,im) = (nanmean(sicEASE(179:182,181,im),1) + ...
%                nanmean(sicEASE(180,180:183,im),2))/2 ;
%            sicEASE(181,181,im) = (nanmean(sicEASE(179:182,181,im),1) + ...
%                nanmean(sicEASE(181,180:183,im),2))/2 ;
%            sicEASE(180,182,im) = (nanmean(sicEASE(179:182,182,im),1) + ...
%                nanmean(sicEASE(180,180:183,im),2))/2 ;
%            sicEASE(181,182,im) = (nanmean(sicEASE(179:182,182,im),1) + ...
%                nanmean(sicEASE(181,180:183,im),2))/2 ;
%            
%            %% Sea-ice velocities
%            % Read in the CESM data
%            uice = double(ncread(filenameU,'uvel',[1 1 mon],[inf inf 1])) ;
%            vice = double(ncread(filenameV,'vvel',[1 1 mon],[inf inf 1])) ;
%            
%            % Rotate the velocity vectors
%            degree = 0 ; % angleU is in radians
%            [uiceR,viceR] = rot_PDR(uice,vice,angleU,degree) ;
%            
%            % Create interpolant function
%            Fu_lat = scatteredInterpolant(latU(:),lonU(:),uiceR(:)) ;
%            Fv_lon = scatteredInterpolant(latU(:),lonU(:),viceR(:)) ;
%            
%            % Interpolate to the EASE Grid
%            uiceEASE(:,:,im) = Fu_lat(LAT, LON) ;
%            viceEASE(:,:,im) = Fv_lon(LAT, LON) ;
%        
%            % Discontinuity problem with uice
%            uiceEASE(181,1:178,im) = NaN ;
%            uiceEASE(180,1:59,im) = NaN ;
%            uiceEASE(182,1:59,im) = NaN ;
%            uiceEASE(180:182,179:181,im) = NaN ;
%            
%            % Fill the line along 180 degree longitude by averaging with
%            % values on both side
%            uiceEASE(180,1:59,im) = nanmean(uiceEASE(178:183,1:59,im),1) ;
%            uiceEASE(182,1:59,im) = nanmean(uiceEASE(180:184,1:59,im),1) ;
%            uiceEASE(181,1:178,im) = nanmean(uiceEASE(180:182,1:178,im),1) ;
%            % Fill points around the North Pole with linear interpolation
%            uiceEASE(160:200,160:200,im) = fillmissing(uiceEASE(160:200,160:200,im),'linear') ;
%            
%            % Discontinuity problem with vice
%            viceEASE(181,1:178,im) = NaN ;
%            viceEASE(180,1:59,im) = NaN ;
%            viceEASE(182,1:59,im) = NaN ;
%            viceEASE(180:182,179:181,im) = NaN ;
%            
%            % Fill the line along 180 degree longitude by averaging with
%            % values on both side
%            viceEASE(180,1:59,im) = nanmean(viceEASE(178:183,1:59,im),1) ;
%            viceEASE(182,1:59,im) = nanmean(viceEASE(180:184,1:59,im),1) ;
%            viceEASE(181,1:178,im) = nanmean(viceEASE(180:182,1:178,im),1) ;
%            % Fill points around the North Pole with linear interpolation
%            viceEASE(160:200,160:200,im) = fillmissing(viceEASE(160:200,160:200,im),'linear') ;
%            
%            % Again, need to rotate the velocity vectors to be able to use
%            % them with the EASE Grid cartesian coordinates
%            degree = 1 ; % Lon is in degrees
%            [uiceEASExy(:,:,im),viceEASExy(:,:,im)] = rot_PDR(uiceEASE(:,:,im),...
%                viceEASE(:,:,im),LON,degree) ;
%            
%            % Change v to positive from top to bottom (+y direction) -> see
%            % comment section at the top
%            viceEASExy(:,:,im) = -viceEASExy(:,:,im) ;
%            
%            % Discontinuity problem at the North Pole again after rotation
%            uiceEASExy(180:182,179:181,im) = NaN ;
%            viceEASExy(180:182,179:181,im) = NaN ;
%            uiceEASExy(160:200,160:200,im) = fillmissing(uiceEASExy(160:200,160:200,im),'linear') ;
%            viceEASExy(160:200,160:200,im) = fillmissing(viceEASExy(160:200,160:200,im),'linear') ;
%            
%            % Declare velocity fields for interpolants
%            uwork = uiceEASExy(:,:,im) ;
%            vwork = viceEASExy(:,:,im) ;
%            
%            % Set all NaNs to zeros to avoid discrepancies with scalar fields
%            uwork(isnan(uwork)) = 0 ;
%            vwork(isnan(vwork)) = 0 ;
%            
%            % Create interpolant function 
%            Fu = scatteredInterpolant(XX(:),YY(:),uwork(:)) ;
%            Fv = scatteredInterpolant(XX(:),YY(:),vwork(:)) ;
%            
%            % Save interpolant
%            save(['IceMotionInterpolant_CESMLE_' member{imem} '_' ...
%                sprintf('%02d',im) '.' int2str(iy) '.mat'],'Fu','Fv') ;
%        end
%        
%        % Set final variables
%        uiceEASE = uiceEASExy ;
%        viceEASE = viceEASExy ;
%        
%        % Save files
%        save(['sicEASE_CESMLE_' member{imem} '_' int2str(iy) '.mat'],...
%            'sicEASE','LAT','LON') ;
%        save(['uiceEASE_CESMLE_' member{imem} '_' int2str(iy) '.mat'],...
%            'uiceEASE','LAT','LON') ;
%        save(['viceEASE_CESMLE_' member{imem} '_' int2str(iy) '.mat'],...
%            'viceEASE','LAT','LON') ;
%    end
%    
%    % Move files
%    movefile(['sicEASE_CESMLE_' member{imem} '_*'],[CESMLEEASEGridMPath ...
%        member{imem} '/']) ;
%    movefile(['uiceEASE_CESMLE_' member{imem} '_*'],[CESMLEEASEGridMPath ...
%        member{imem} '/']) ;
%    movefile(['viceEASE_CESMLE_' member{imem} '_*'],[CESMLEEASEGridMPath ...
%        member{imem} '/']) ;
%    movefile(['IceMotionInterpolant_CESMLE_' member{imem} '_*'],...
%        [CESMLEInterpMPath member{imem} '/']) ;
%    toc

    %% 10-m wind
%     % Read in the CESM data
%     filenameU10 = ([CESMLERawMPath 'RCP8.5/U10/b.e11.BRCP85C5CNBDRD.f09_g16.' ...
%         member{imem} '.cam.h0.U10.200601-210012.nc']) ;
%     ystart = 2006 ;
% 
%     latU10 = double(ncread(filenameU10,'lat',97,inf)) ;
%     lonU10 = double(ncread(filenameU10,'lon',1,inf)) ;
%     lonU10(lonU10>180) = lonU10(lonU10>180)-360 ;
% 
%     [LATU10,LONU10] = meshgrid(latU10,lonU10) ;
% 
%     % Declare initial variable
%     U10EASE = NaN(361,361,12) ;
% 
%     % Looping through the years
%     for iy = 2006:2100
% 
%         % Looping through the months
%         for im = 1:12
%             disp(['Year: ' int2str(iy) ', Month: ' int2str(im)])
%             mon = (iy-ystart)*12 + im ;
% 
%             % Read in the CESM data
%             U10 = double(ncread(filenameU10,'U10',[1 97 mon],[inf inf 1])) ;
% 
% %             figure(1), clf
% %             pcolor(U10'); shading flat; colorbar
% 
%             % Create interpolant function
%             FU10 = scatteredInterpolant(LATU10(:),LONU10(:),U10(:)) ;
% 
%             % Interpolate to the EASE Grid
%             U10EASE(:,:,im) = FU10(LAT,LON) ;
% 
% %             figure(2), clf
% %             pcolor(XX,YY,U10EASE(:,:,im)); shading flat; colorbar; axis ij
% 
%             % Discontinuity problem
%             U10EASE(180,1:180,im) = NaN ;
%             U10EASE(180,1:180,im) = nanmean(U10EASE(179:181,1:180,im),1) ;
% 
% %             figure(3), clf
% %             pcolor(XX,YY,U10EASE(:,:,im)); shading flat; colorbar; axis ij
% %             pause
%         end
% 
%         % Save file
%         save(['U10EASE_CESMLE_' member{imem} '_' int2str(iy) '.mat'],...
%             'U10EASE','LAT','LON') ;
%     end
% 
%     % Move file
%     movefile(['U10EASE_CESMLE_' member{imem} '_*'],[CESMLEEASEGridMPath ...
%         member{imem} '/']) ;
%     toc
    
    %% Ice thickness categories
%     
%         %% Category 1 (0-0.60m)
%         % Read in the CESM data
%         filenameA = ([CESMLERawMPath 'RCP8.5/aice001/b.e11.BRCP85C5CNBDRD.f09_g16.' ...
%             member{imem} '.cice.h.aicen001_nh.200601-210012.nc']) ;
%         
%         % Declare initial variables
%         aice1EASE = NaN(361,361,12) ;
%         ystart = 2006 ;
% 
%         % Looping through the years
%         for iy = 2006:2100
%             
%             % Looping through the months
%             for im = 1:12
%                 disp(['Year: ' int2str(iy) ', Month: ' int2str(im)])
%                 mon = (iy-ystart)*12 + im ;
%                 
%                 % Read in the CESM data
%                 aicen = ncread(filenameA,'aicen001',[1 1 mon],[inf inf 1]) ;
%                 
%                 % Create interpolant function
%                 Fa = scatteredInterpolant(latT(:), lonT(:), aicen(:)) ;
%                 
%                 % Interpolate to the EASE Grid
%                 aice1EASE(:,:,im) = Fa(LAT,LON) ;
%                 
%                 % Discontinuity problem
%                 aice1EASE(181,1:180,im) = nanmean(aice1EASE(180:182,1:180,im),1) ;
%             end
%             
%             % Save file
%             save(['aice1EASE_' int2str(iy) '.mat'], 'aice1EASE', 'LAT', 'LON') ;
%         end
%         
%         %% Category 2 (0.60-1.40m)
%         % Read in the CESM data
%         filenameA = ([CESMLERawMPath 'RCP8.5/aice002/b.e11.BRCP85C5CNBDRD.f09_g16.' ...
%             member{imem} '.cice.h.aicen002_nh.200601-210012.nc']) ;
%         
%         % Declare initial variables
%         aice2EASE = NaN(361,361,12) ;
%         ystart = 2006 ;
% 
%         % Looping through the years
%         for iy = 2006:2100
%             
%             % Looping through the months
%             for im = 1:12
%                 disp(['Year: ' int2str(iy) ', Month: ' int2str(im)])
%                 mon = (iy-ystart)*12 + im ;
%                 
%                 % Read in the CESM data
%                 aicen = ncread(filenameA,'aicen002',[1 1 mon],[inf inf 1]) ;
%                 
%                 % Create interpolant function
%                 Fa = scatteredInterpolant(latT(:), lonT(:), aicen(:)) ;
%                 
%                 % Interpolate to the EASE Grid
%                 aice2EASE(:,:,im) = Fa(LAT,LON) ;
%                 
%                 % Discontinuity problem
%                 aice2EASE(181,1:180,im) = nanmean(aice2EASE(180:182,1:180,im),1) ;
%             end
%             
%             % Save file
%             save(['aice2EASE_' int2str(iy) '.mat'], 'aice2EASE', 'LAT', 'LON') ;
%         end
%         
%         %% Category 3 (1.40-2.40m)
%         % Read in the CESM data
%         filenameA = ([CESMLERawMPath 'RCP8.5/aice003/b.e11.BRCP85C5CNBDRD.f09_g16.' ...
%             member{imem} '.cice.h.aicen003_nh.200601-210012.nc']) ;
%         
%         % Declare initial variables
%         aice3EASE = NaN(361,361,12) ;
%         ystart = 2006 ;
% 
%         % Looping through the years
%         for iy = 2006:2100
%             
%             % Looping through the months
%             for im = 1:12
%                 disp(['Year: ' int2str(iy) ', Month: ' int2str(im)])
%                 mon = (iy-ystart)*12 + im ;
%                 
%                 % Read in the CESM data
%                 aicen = ncread(filenameA,'aicen003',[1 1 mon],[inf inf 1]) ;
%                 
%                 % Create interpolant function
%                 Fa = scatteredInterpolant(latT(:), lonT(:), aicen(:)) ;
%                 
%                 % Interpolate to the EASE Grid
%                 aice3EASE(:,:,im) = Fa(LAT,LON) ;
%                 
%                 % Discontinuity problem
%                 aice3EASE(181,1:180,im) = nanmean(aice3EASE(180:182,1:180,im),1) ;
%             end
%             
%             % Save file
%             save(['aice3EASE_' int2str(iy) '.mat'], 'aice3EASE', 'LAT', 'LON') ;
%         end
%         
%         %% Category 4 (2.40-3.60m)
%         % Read in the CESM data
%         filenameA = ([CESMLERawMPath 'RCP8.5/aice004/b.e11.BRCP85C5CNBDRD.f09_g16.' ...
%             member{imem} '.cice.h.aicen004_nh.200601-210012.nc']) ;
%         
%         % Declare initial variables
%         aice4EASE = NaN(361,361,12) ;
%         ystart = 2006 ;
% 
%         % Looping through the years
%         for iy = 2006:2100
%             
%             % Looping through the months
%             for im = 1:12
%                 disp(['Year: ' int2str(iy) ', Month: ' int2str(im)])
%                 mon = (iy-ystart)*12 + im ;
%                 
%                 % Read in the CESM data
%                 aicen = ncread(filenameA,'aicen004',[1 1 mon],[inf inf 1]) ;
%                 
%                 % Create interpolant function
%                 Fa = scatteredInterpolant(latT(:), lonT(:), aicen(:)) ;
%                 
%                 % Interpolate to the EASE Grid
%                 aice4EASE(:,:,im) = Fa(LAT,LON) ;
%                 
%                 % Discontinuity problem
%                 aice4EASE(181,1:180,im) = nanmean(aice4EASE(180:182,1:180,im),1) ;
%             end
%             
%             % Save file
%             save(['aice4EASE_' int2str(iy) '.mat'], 'aice4EASE', 'LAT', 'LON') ;
%         end
%         
%         %% Category 5 (>3.60m)
%         % Read in the CESM data
%         filenameA = ([CESMLERawMPath 'RCP8.5/aice005/b.e11.BRCP85C5CNBDRD.f09_g16.' ...
%             member{imem} '.cice.h.aicen005_nh.200601-210012.nc']) ;
%         
%         % Declare initial variables
%         aice5EASE = NaN(361,361,12) ;
%         ystart = 2006 ;
% 
%         % Looping through the years
%         for iy = 2006:2100
%             
%             % Looping through the months
%             for im = 1:12
%                 disp(['Year: ' int2str(iy) ', Month: ' int2str(im)])
%                 mon = (iy-ystart)*12 + im ;
%                 
%                 % Read in the CESM data
%                 aicen = ncread(filenameA,'aicen005',[1 1 mon],[inf inf 1]) ;
%                 
%                 % Create interpolant function
%                 Fa = scatteredInterpolant(latT(:), lonT(:), aicen(:)) ;
%                 
%                 % Interpolate to the EASE Grid
%                 aice5EASE(:,:,im) = Fa(LAT,LON) ;
%                 
%                 % Discontinuity problem
%                 aice5EASE(181,1:180,im) = nanmean(aice5EASE(180:182,1:180,im),1) ;
%             end
%             
%             % Save file
%             save(['aice5EASE_' int2str(iy) '.mat'], 'aice5EASE', 'LAT', 'LON') ;
%         end
%     
%     % Move files
%     movefile('aice*',[CESMLEEASEGridMPath member{imem} '/']) ;
%     toc
    
    %% Sea-ice thickness
%     % Read in the CESM data
%     filenameT = ([CESMLERawMPath 'RCP8.5/SIT/b.e11.BRCP85C5CNBDRD.f09_g16.' ...
%         member{imem} '.cice.h.hi_nh.200601-210012.nc']) ;
%     
%     % Declare initial variables
%     sitEASE = NaN(361,361,12) ;
%     ystart = 2006 ;
%     
%     % Looping through the years
%     for iy = 2006:2100
%         
%         % Looping through the months
%         for im = 1:12
%             disp(['Year: ' int2str(iy) ', Month: ' int2str(im)])
%             mon = (iy-ystart)*12 + im ;
%             
%             % Read in the CESM data
%             sit = ncread(filenameT,'hi',[1 1 mon],[inf inf 1]) ;
%             
%             % Create interpolant function
%             Ft = scatteredInterpolant(latT(:), lonT(:), sit(:)) ;
%             
%             % Interpolate to the EASE Grid
%             sitEASE(:,:,im) = Ft(LAT,LON) ;
%             
%             % Discontinuity problem
%             sitEASE(181,1:180,im) = nanmean(sitEASE(180:182,1:180,im),1) ;
%         end
%         
%         % Save file
%         save(['sitEASE_' int2str(iy) '.mat'], 'sitEASE', 'LAT', 'LON') ;
%     end
%     
%     % Move file
%     movefile('sitEASE*',[CESMLEEASEGridMPath member{imem} '/']) ;
%     toc

    %% Sea-level Pressure
%     % Read in the CESM data
%     filenameP = ([CESMLERawMPath 'RCP8.5/PSL/b.e11.BRCP85C5CNBDRD.f09_g16.' ...
%         member{imem} '.cam.h0.PSL.200601-210012.nc']) ;
% 
%     % CESM-grid lats and lons
%     % for CESM PSL files: Lat = 129-192 = north of 30N (matches 361x361 grid)
%     latP_tmp = double(ncread(filenameP,'lat',129,inf)) ;
%     lonP_tmp = double(ncread(filenameP,'lon',1,inf)) ;
%     lonP_tmp(lonP_tmp>180) = lonP_tmp(lonP_tmp>180)-360 ; % change to lon from -180 to 180
% 
%     % Create lat/lon matrices
%     latP = repmat(latP_tmp',max(size(lonP_tmp)),1) ;
%     lonP = repmat(lonP_tmp,1,max(size(latP_tmp))) ;
%     
%     % Declare initial variables
%     slpEASE = NaN(361,361,12) ;
%     ystart = 2006 ;
%     
%     % Looping through the years
%     for iy = 2006:2100
%         
%         % Looping through the months
%         for im = 1:12
%             disp(['Year: ' int2str(iy) ', Month: ' int2str(im)])
%             mon = (iy-ystart)*12 + im ;
%             
%             % Read in the CESM data
%             SLP = double(ncread(filenameP,'PSL',[1 129 mon],[inf inf 1])) ; % in Pa
%             
%             % Create interpolant function
%             Fslp = scatteredInterpolant(latP(:),lonP(:),SLP(:)) ;
%             
%             % Interpolate to the EASE Grid
%             slpEASE(:,:,im) = Fslp(LAT, LON) ;
%             slpEASE(:,:,im) = slpEASE(:,:,im)./100 ; % Put SLP in hPa
%             
% %             figure(1), clf
% %             pcolor(XX,YY,slpEASE(:,:,im)); shading flat
% %             hold on, box on, axis ij equal tight
% %             colorbar
% %             caxis([995 1025])
% %             S = surf(XX,YY,new_CESMmask) ; shading flat
% %             z = get(S,'ZData') ;
% %             set(S,'FaceColor',[0.8 0.8 0.8],'ZData',z+10)
% %             xlim([181-l 181+l])
% %             ylim([181-l 181+l])
% % %             set(gca,'XTick',[],'YTick',[])
% %             pause
%             
% %             % Discontinuity problem
% %             slpEASE(178,1:43,im) = nanmean(slpEASE(177:181,1:43,im),1) ;
% %             slpEASE(179,1:89,im) = nanmean(slpEASE(178:181,1:89,im),1) ;
% %             slpEASE(180,1:135,im) = nanmean(slpEASE(179:181,1:135,im),1) ;
% %             
% %             figure(2), clf
% %             pcolor(XX,YY,slpEASE(:,:,im)); shading flat
% %             hold on, box on, axis ij equal tight
% %             colorbar
% %             caxis([995 1025])
% %             S = surf(XX,YY,new_CESMmask) ; shading flat
% %             z = get(S,'ZData') ;
% %             set(S,'FaceColor',[0.8 0.8 0.8],'ZData',z+10)
% %             xlim([181-l 181+l])
% %             ylim([181-l 181+l])
% %             set(gca,'XTick',[],'YTick',[])
% %             pause
%         end
%         
%         % Save file
%         save(['slpEASE_CESMLE_' member{imem} '_' int2str(iy) '.mat'], 'slpEASE', 'LAT', 'LON') ;
%     end
%     
%     % Move file
%     movefile(['slpEASE_CESMLE_' member{imem} '_*.mat'],...
%         [CESMLEEASEGridMPath member{imem} '/']) ;
%     toc
    
    %% 2-m air temperature
%     % Read in the CESM data
%     filenameTEMP = ([CESMLERawMPath 'RCP8.5/TREFHT/b.e11.BRCP85C5CNBDRD.f09_g16.' ...
%         member{imem} '.cam.h0.TREFHT.200601-210012.nc']) ;
% 
%     % CESM-grid lats and lons
%     % for CESM TREFHT files: Lat = 129-192 = north of 30N (matches 361x361 grid)
%     latTEMP_tmp = double(ncread(filenameTEMP,'lat',129,inf)) ;
%     lonTEMP_tmp = double(ncread(filenameTEMP,'lon',1,inf)) ;
%     lonTEMP_tmp(lonTEMP_tmp>180) = lonTEMP_tmp(lonTEMP_tmp>180)-360 ; % change to lon from -180 to 180
% 
%     % Create lat/lon matrices
%     latTEMP = repmat(latTEMP_tmp',max(size(lonTEMP_tmp)),1) ;
%     lonTEMP = repmat(lonTEMP_tmp,1,max(size(latTEMP_tmp))) ;
%     
%     % Declare initial variables
%     tempEASE = NaN(361,361,12) ;
%     ystart = 2006 ;
%     
%     % Looping through the years
%     for iy = 2006:2100
%         
%         % Looping through the months
%         for im = 1:12
%             disp(['Year: ' int2str(iy) ', Month: ' int2str(im)])
%             mon = (iy-ystart)*12 + im ;
%             
%             % Read in the CESM data
%             temp = double(ncread(filenameTEMP,'TREFHT',[1 129 mon],[inf inf 1])) ; % in K
%             
%             % Create interpolant function
%             Ftemp = scatteredInterpolant(latTEMP(:),lonTEMP(:),temp(:)) ;
%             
%             % Interpolate to the EASE Grid
%             tempEASE(:,:,im) = Ftemp(LAT, LON) ;
%             
% %             figure(1), clf
% %             pcolor(XX,YY,tempEASE(:,:,im)); shading flat
% %             hold on, box on, axis ij equal tight
% %             colorbar
% %             caxis([220 300])
% %             S = surf(XX,YY,new_CESMmask) ; shading flat
% %             z = get(S,'ZData') ;
% %             set(S,'FaceColor',[0.8 0.8 0.8],'ZData',z+10)
% %             xlim([181-l 181+l])
% %             ylim([181-l 181+l])
% % %             set(gca,'XTick',[],'YTick',[])
% %             pause
%             
% %             % Discontinuity problem
% %             tempEASE(178,1:43,im) = nanmean(tempEASE(177:181,1:43,im),1) ;
% %             tempEASE(179,1:89,im) = nanmean(tempEASE(178:181,1:89,im),1) ;
% %             tempEASE(180,1:135,im) = nanmean(tempEASE(179:181,1:135,im),1) ;
% %             
% %             figure(2), clf
% %             pcolor(XX,YY,tempEASE(:,:,im)); shading flat
% %             hold on, box on, axis ij equal tight
% %             colorbar
% %             caxis([995 1025])
% %             S = surf(XX,YY,new_CESMmask) ; shading flat
% %             z = get(S,'ZData') ;
% %             set(S,'FaceColor',[0.8 0.8 0.8],'ZData',z+10)
% %             xlim([181-l 181+l])
% %             ylim([181-l 181+l])
% %             set(gca,'XTick',[],'YTick',[])
% %             pause
%         end
%         
%         % Save file
%         save(['tempEASE_CESMLE_' member{imem} '_' int2str(iy) '.mat'], 'tempEASE', 'LAT', 'LON') ;
%     end
%     
%     % Move file
%     movefile(['tempEASE_CESMLE_' member{imem} '_*.mat'],...
%         [CESMLEEASEGridMPath member{imem} '/']) ;
%     toc

    %% U and V component of the wind
    % Read in the CESM data
    filenameUWIND = ([CESMLERawMPath 'RCP8.5/U/b.e11.BRCP85C5CNBDRD.f09_g16.' ...
       member{imem} '.cam.h0.U.200601-210012.nc']) ; 
    filenameVWIND = ([CESMLERawMPath 'RCP8.5/V/b.e11.BRCP85C5CNBDRD.f09_g16.' ...
       member{imem} '.cam.h0.V.200601-210012.nc']) ;
   
    % CESM-grid lats and lons
    % for CESM U and V files: Lat = 129-192 = north of 30N (matches 361x361 grid)
    latWIND_tmp = double(ncread(filenameUWIND,'lat',129,inf)) ;
    lonWIND_tmp = double(ncread(filenameUWIND,'lon',1,inf)) ;
    lonWIND_tmp(lonWIND_tmp>180) = lonWIND_tmp(lonWIND_tmp>180)-360 ; % change to lon from -180 to 180

    % Create lat/lon matrices
    latWIND = repmat(latWIND_tmp',max(size(lonWIND_tmp)),1) ;
    lonWIND = repmat(lonWIND_tmp,1,max(size(latWIND_tmp))) ;
   
    % Declare initial variables
    uwindEASE   = NaN(361,361,12) ;
    vwindEASE   = NaN(361,361,12) ;
    uwindEASExy = NaN(361,361,12) ;
    vwindEASExy = NaN(361,361,12) ;
    ystart = 2006 ; 
   
    for iy = 2006:2100
        disp(['Year: ' int2str(iy)])
       
        % Looping through the months
        for im = 1:12
            disp(['Year: ' int2str(iy) ', Month: ' int2str(im)])
            mon = (iy-ystart)*12 + im ;

            % Read in the CESM data
            % Lat = 129-192 = north of 30N (matches 361x361 grid)
            % Last level (30) = surface (defined positive down, starts at TOA)
            uwind = double(ncread(filenameUWIND,'U',[1 129 30 mon],[inf inf 1 1])) ; % in m/s
            vwind = double(ncread(filenameVWIND,'V',[1 129 30 mon],[inf inf 1 1])) ; % in m/s

            % Create interpolant function
            Fuwind_lat = scatteredInterpolant(latWIND(:),lonWIND(:),uwind(:)) ;
            Fvwind_lon = scatteredInterpolant(latWIND(:),lonWIND(:),vwind(:)) ;

            % Interpolate to the EASE Grid
            uwindEASE(:,:,im) = Fuwind_lat(LAT,LON) ;
            vwindEASE(:,:,im) = Fvwind_lon(LAT,LON) ;
            
%             figure(1), clf
%             contour(XX,YY,bathyEASE,[0,0],'LineColor','k','LineWidth',2)
%             hold on, box on, axis ij equal tight
%             quiver(XX(1:5:361,1:5:361),YY(1:5:361,1:5:361),uwindEASE(1:5:361,1:5:361,im),vwindEASE(1:5:361,1:5:361,im),2)
%             P = pcolor(XX,YY,vwindEASE(:,:,im)) ; shading flat
%             z = get(P,'ZData') ;
%             set(P,'ZData',z-10)
%             colorbar

            % Again, need to rotate the velocity vectors to be able to use
            % them with the EASE Grid cartesian coordinates
            degree = 1 ; % Lon is in degrees
            [uwindEASExy(:,:,im),vwindEASExy(:,:,im)] = rot_PDR(uwindEASE(:,:,im),...
               vwindEASE(:,:,im),LON,degree) ;
            
            % Discontinuity problem at the North Pole after rotation
            uwindEASExy(181,181,im) = NaN ;
            vwindEASExy(181,181,im) = NaN ;

            % Fill point at the North Pole with linear interpolation
            uwindEASExy(180:182,180:182,im) = fillmissing(uwindEASExy(180:182,180:182,im),'linear') ;
            vwindEASExy(180:182,180:182,im) = fillmissing(vwindEASExy(180:182,180:182,im),'linear') ;
            
%             figure(2), clf
%             contour(XX,YY,bathyEASE,[0,0],'LineColor','k','LineWidth',2)
%             hold on, box on, axis ij equal tight
%             quiver(XX(1:5:361,1:5:361),YY(1:5:361,1:5:361),uwindEASExy(1:5:361,1:5:361,im),vwindEASExy(1:5:361,1:5:361,im),2)
%             P = pcolor(XX,YY,vwindEASExy(:,:,im)) ; shading flat
%             z = get(P,'ZData') ;
%             set(P,'ZData',z-10)
%             colorbar
%             pause

            % Change v to positive from top to bottom (+y direction) -> see
            % comment section at the top
            vwindEASExy(:,:,im) = -vwindEASExy(:,:,im) ;
            
%             figure(3), clf
%             contour(XX,YY,bathyEASE,[0,0],'LineColor','k','LineWidth',2)
%             hold on, box on, axis ij equal tight
%             quiver(XX(1:5:361,1:5:361),YY(1:5:361,1:5:361),uwindEASExy(1:5:361,1:5:361,im),vwindEASExy(1:5:361,1:5:361,im),2)
%             P = pcolor(XX,YY,vwindEASExy(:,:,im)) ; shading flat
%             z = get(P,'ZData') ;
%             set(P,'ZData',z-10)
%             colorbar
%             pause

        end
       
        % Set final variables
        uwindEASE = uwindEASExy ;
        vwindEASE = vwindEASExy ;
       
        % Save files
        save(['uwindEASE_CESMLE_' member{imem} '_' int2str(iy) '.mat'],...
           'uwindEASE','LAT','LON') ;
        save(['vwindEASE_CESMLE_' member{imem} '_' int2str(iy) '.mat'],...
           'vwindEASE','LAT','LON') ;
    end

    % Move files
    movefile(['uwindEASE_CESMLE_' member{imem} '_*'],[CESMLEEASEGridMPath ...
       member{imem} '/']) ;
    movefile(['vwindEASE_CESMLE_' member{imem} '_*'],[CESMLEEASEGridMPath ...
       member{imem} '/']) ;
    toc

end
