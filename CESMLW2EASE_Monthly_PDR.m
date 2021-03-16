% CESMLW2EASE_Monthly_PDR
%
% A script to:
%   * Load CESM 2 degree C model output
%   * Interpolate variables (aice, uvel, vvel, PSL) onto EASE Grid
%   * Rotate the velocity vectors towards the real North Pole
%   * Interpolate onto EASE Grid
%   * Store the interpolants and fields in a set of monthly files.
%
% Patricia, May 2018
% Last update: May 2020

%% Common variables
% % When working on Mac
% FAITPaths_Mac_PDR

% When working on Cheyenne
FAITPaths_Cheyenne

% Ensemble members
member = {'001' '002' '003' '004' '005' '006' '007' '008' '009' '010' '011'} ;

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
clear north_x_y_lat_lon

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
    
%     % Read in the CESM data
%     filenameC = ([CESMLWRawMPath 'aice/b.e11.BRCP26C5CNBDRD.f09_g16.2pt0degC.' ...
%         member{imem} '.cice.h.aice_nh.200601-210012.nc']) ;
%     filenameU = ([CESMLWRawMPath 'uvel/b.e11.BRCP26C5CNBDRD.f09_g16.2pt0degC.' ...
%         member{imem} '.cice.h.uvel_nh.200601-210012.nc']) ;
%     filenameV = ([CESMLWRawMPath 'vvel/b.e11.BRCP26C5CNBDRD.f09_g16.2pt0degC.' ...
%         member{imem} '.cice.h.vvel_nh.200601-210012.nc']) ;
%     
%     % Declare initial variables
%     sicEASE = NaN(361,361,12) ;
%     uiceEASE = NaN(361,361,12) ;
%     viceEASE = NaN(361,361,12) ;
%     uiceEASExy = NaN(361,361,12) ;
%     viceEASExy = NaN(361,361,12) ;
%     ystart = 2006 ; 
%     
%     %% Looping through the years
%     for iy = 2006:2100
%         disp(['Year: ' int2str(iy)])
%         
%         % Looping through the months
%         for im = 1:12
% %             disp(['Year: ' int2str(iy) ', Month: ' int2str(im)])
%             mon = (iy-ystart)*12 + im ;
%             
%             %% Sea-ice concentration
%             % Read in the CESM data
%             sic = ncread(filenameC,'aice',[1 1 mon],[inf inf 1]) ;
%             
%             % Create interpolant function
%             Fc = scatteredInterpolant(latT(:), lonT(:), sic(:)) ;
%             
%             % Interpolate to the EASE Grid
%             sicEASE(:,:,im) = Fc(LAT,LON) ;
%             
%             % Discontinuity problem
%             sicEASE(181,1:180,im) = NaN ;
%             sicEASE(180,1:59,im) = NaN ;
%             sicEASE(182,1:59,im) = NaN ;
%             sicEASE(180,180,im) = NaN ;
%             sicEASE(182,181,im) = NaN ;
%             sicEASE(180,181,im) = NaN ;
%             sicEASE(181,181,im) = NaN ;
%             sicEASE(180,182,im) = NaN ;
%             sicEASE(181,182,im) = NaN ;
%             
%             sicEASE(180,1:59,im) = nanmean(sicEASE(178:183,1:59,im),1) ;
%             sicEASE(182,1:59,im) = nanmean(sicEASE(180:184,1:59,im),1) ;
%             sicEASE(181,1:180,im) = nanmean(sicEASE(180:182,1:180,im),1) ;
%             sicEASE(180,180,im) = nanmean(sicEASE(179:181,180,im),1) ;
%             sicEASE(182,181,im) = nanmean(sicEASE(182,180:182,im),2) ;
%             sicEASE(180,181,im) = (nanmean(sicEASE(179:182,181,im),1) + ...
%                 nanmean(sicEASE(180,180:183,im),2))/2 ;
%             sicEASE(181,181,im) = (nanmean(sicEASE(179:182,181,im),1) + ...
%                 nanmean(sicEASE(181,180:183,im),2))/2 ;
%             sicEASE(180,182,im) = (nanmean(sicEASE(179:182,182,im),1) + ...
%                 nanmean(sicEASE(180,180:183,im),2))/2 ;
%             sicEASE(181,182,im) = (nanmean(sicEASE(179:182,182,im),1) + ...
%                 nanmean(sicEASE(181,180:183,im),2))/2 ;
%             
%             %% Sea-ice velocities
%             % Read in the CESM data
%             uice = double(ncread(filenameU,'uvel',[1 1 mon],[inf inf 1])) ;
%             vice = double(ncread(filenameV,'vvel',[1 1 mon],[inf inf 1])) ;
%             
%             % Rotate the velocity vectors
%             degree = 0 ; % angleU is in radians
%             [uiceR,viceR] = rot_PDR(uice,vice,angleU,degree) ;
%             
%             % Create interpolant function
%             Fu_lat = scatteredInterpolant(latU(:),lonU(:),uiceR(:)) ;
%             Fv_lon = scatteredInterpolant(latU(:),lonU(:),viceR(:)) ;
%             
%             % Interpolate to the EASE Grid
%             uiceEASE(:,:,im) = Fu_lat(LAT, LON) ;
%             viceEASE(:,:,im) = Fv_lon(LAT, LON) ;
%             
%             % Discontinuity problem with uice
%             uiceEASE(181,1:178,im) = NaN ;
%             uiceEASE(180,1:59,im) = NaN ;
%             uiceEASE(182,1:59,im) = NaN ;
%             uiceEASE(180:182,179:181,im) = NaN ;
%             
%             % Fill the line along 180 degree longitude by averaging with
%             % values on both side
%             uiceEASE(180,1:59,im) = nanmean(uiceEASE(178:183,1:59,im),1) ;
%             uiceEASE(182,1:59,im) = nanmean(uiceEASE(180:184,1:59,im),1) ;
%             uiceEASE(181,1:178,im) = nanmean(uiceEASE(180:182,1:178,im),1) ;
%             % Fill points around the North Pole with linear interpolation
%             uiceEASE(160:200,160:200,im) = fillmissing(uiceEASE(160:200,160:200,im),'linear') ;
%                         
%             % Discontinuity problem with vice
%             viceEASE(181,1:178,im) = NaN ;
%             viceEASE(180,1:59,im) = NaN ;
%             viceEASE(182,1:59,im) = NaN ;
%             viceEASE(180:182,179:181,im) = NaN ;
%             
%             % Fill the line along 180 degree longitude by averaging with
%             % values on both side
%             viceEASE(180,1:59,im) = nanmean(viceEASE(178:183,1:59,im),1) ;
%             viceEASE(182,1:59,im) = nanmean(viceEASE(180:184,1:59,im),1) ;
%             viceEASE(181,1:178,im) = nanmean(viceEASE(180:182,1:178,im),1) ;
%             % Fill points around the North Pole with linear interpolation
%             viceEASE(160:200,160:200,im) = fillmissing(viceEASE(160:200,160:200,im),'linear') ;
%             
%             % Again, need to rotate the velocity vectors to be able to use
%             % them with the EASE Grid cartesian coordinates
%             degree = 1 ; % Lon is in degrees
%             [uiceEASExy(:,:,im),viceEASExy(:,:,im)] = rot_PDR(uiceEASE(:,:,im),...
%                 viceEASE(:,:,im),LON,degree) ;
%             
%             % Change v to positive from top to bottom (+y direction) -> see
%             % comment section at the top
%             viceEASExy(:,:,im) = -viceEASExy(:,:,im) ;
%             
%             % Discontinuity problem at the North Pole again after rotation
%             uiceEASExy(180:182,179:181,im) = NaN ;
%             viceEASExy(180:182,179:181,im) = NaN ;
%             uiceEASExy(160:200,160:200,im) = fillmissing(uiceEASExy(160:200,160:200,im),'linear') ;
%             viceEASExy(160:200,160:200,im) = fillmissing(viceEASExy(160:200,160:200,im),'linear') ;
%             
%             % Declare velocity fields for interpolants
%             uwork = uiceEASExy(:,:,im) ;
%             vwork = viceEASExy(:,:,im) ;
%             
%             % Set all NaNs to zeros to avoid discrepancies with scalar fields
%             uwork(isnan(uwork)) = 0 ;
%             vwork(isnan(vwork)) = 0 ;
%             
%             % Create interpolant function 
%             Fu = scatteredInterpolant(XX(:),YY(:),uwork(:)) ;
%             Fv = scatteredInterpolant(XX(:),YY(:),vwork(:)) ;
%             
%             % Save interpolant
%             save(['IceMotionInterpolant_CESMLW_' member{imem} '_' ...
%                 sprintf('%02d',im) '.' int2str(iy) '.mat'],'Fu','Fv') ;
%         end
%         
%         % Set final variables
%         uiceEASE = uiceEASExy ;
%         viceEASE = viceEASExy ;
%         
%         % Save files
%         save(['sicEASE_CESMLW_' member{imem} '_' int2str(iy) '.mat'],...
%             'sicEASE','LAT','LON') ;
%         save(['uiceEASE_CESMLW_' member{imem} '_' int2str(iy) '.mat'],...
%         'uiceEASE','LAT','LON') ;
%         save(['viceEASE_CESMLW_' member{imem} '_' int2str(iy) '.mat'],...
%             'viceEASE','LAT','LON') ;
%     end
%     
%     % Move files
%     movefile(['sicEASE_CESMLW_' member{imem} '_*'],[CESMLWEASEGridMPath ...
%         member{imem} '/']) ;
%     movefile(['uiceEASE_CESMLW_' member{imem} '_*'],[CESMLWEASEGridMPath ...
%         member{imem} '/']) ;
%     movefile(['viceEASE_CESMLW_' member{imem} '_*'],[CESMLWEASEGridMPath ...
%         member{imem} '/']) ;
%     movefile(['IceMotionInterpolant_CESMLW_' member{imem} '_*'],...
%         [CESMLWInterpMPath member{imem} '/']) ;
%     toc

    %% Sea-level Pressure
%     % Read in the CESM data
%     filenameP = ([CESMLWRawMPath 'PSL/b.e11.BRCP26C5CNBDRD.f09_g16.2pt0degC.' ...
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
% %             % Discontinuity problem
% %             slpEASE(178,1:43,im) = nanmean(slpEASE(177:181,1:43,im),1) ;
% %             slpEASE(179,1:89,im) = nanmean(slpEASE(178:181,1:89,im),1) ;
% %             slpEASE(180,1:135,im) = nanmean(slpEASE(179:181,1:135,im),1) ;
%         end
%         
%         % Save files
%         save(['slpEASE_CESMLW_' member{imem} '_' int2str(iy) '.mat'],...
%             'slpEASE','LAT','LON') ;
%     end
%     
%     % Move files
%     movefile(['slpEASE_CESMLW_' member{imem} '_*.mat'],[CESMLWEASEGridMPath ...
%         member{imem} '/']) ;
%     toc
    
    %% 2-m air temperature
%     % Read in the CESM data
%     filenameTEMP = ([CESMLWRawMPath 'TREFHT/b.e11.BRCP26C5CNBDRD.f09_g16.2pt0degC.' ...
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
%         save(['tempEASE_CESMLW_' member{imem} '_' int2str(iy) '.mat'], 'tempEASE', 'LAT', 'LON') ;
%     end
%     
%     % Move file
%     movefile(['tempEASE_CESMLW_' member{imem} '_*.mat'],...
%         [CESMLWEASEGridMPath member{imem} '/']) ;
%     toc
    
    %% U and V component of the wind
    % Read in the CESM data
    filenameUWIND = ([CESMLWRawMPath 'U/b.e11.BRCP26C5CNBDRD.f09_g16.2pt0degC.' ...
       member{imem} '.cam.h0.U.200601-210012.nc']) ; 
    filenameVWIND = ([CESMLWRawMPath 'V/b.e11.BRCP26C5CNBDRD.f09_g16.2pt0degC.' ...
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
        save(['uwindEASE_CESMLW_' member{imem} '_' int2str(iy) '.mat'],...
           'uwindEASE','LAT','LON') ;
        save(['vwindEASE_CESMLW_' member{imem} '_' int2str(iy) '.mat'],...
           'vwindEASE','LAT','LON') ;
    end

    % Move files
    movefile(['uwindEASE_CESMLW_' member{imem} '_*'],[CESMLWEASEGridMPath ...
       member{imem} '/']) ;
    movefile(['vwindEASE_CESMLW_' member{imem} '_*'],[CESMLWEASEGridMPath ...
       member{imem} '/']) ;
    toc
    
end