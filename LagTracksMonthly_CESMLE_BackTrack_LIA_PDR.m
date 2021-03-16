%  LagTracksMonthly_CESMLE_BackTrack_LIA_PDR
%
%  A scrip to:
%    * Track all points in the LIA forward or backward in time for a range
%    of startyears and startmonths
%    * Done for all members of the CESM-LE
%    * Save all positions along the tracks in a .mat file
%
%  This script is a modified version of LagTracksEEZMonthly_CESMLE_PDR.m
%  that was used to track ice parcels forward in time for the EEZFuture
%  project (Earth's Future (2020) paper).
%
%  The starting points are created by TromsoPatches_PDR.m and are stored
%  under the directory XYIndices.
%
%  The routine checks that the active points are still in the ice 
%  pack at each time step.    
%
%  Reconfigured the run both forward and backward.  
%
%  Mismatch between the EASE-25 and CESM grids led to ice being 
%  over land.  The scripts (1) used nearest neighbor for checking 
%  ice concentration and (2) pulls ice on land back along the trajectory
%  of its last step to find the nearest ocean point.  
%
%  Patricia DeRepentigny, Robert Newton, Bruno Tremblay, 2020
%
%  Last modified: May 2020

%% Common Variables
% Run external script to set path statements for specific environment.
FAITPaths_IceShedMS

% % Read EASE Grid coordinates
% Depending on whether the EASE grid is already loaded:
% load('north_x_y_lat_lon')
% x       = north_x_y_lat_lon(:,1) + 1 ;
% y       = north_x_y_lat_lon(:,2) + 1 ;
% XX      = reshape(x,361,361) ;
% YY      = reshape(y,361,361) ;
% clear north_x_y_lat_lon

% Ensemble members
member = {'001' '002' '003' '004' '005' '006' '007' '008' '009' '010' ...
    '011' '012' '013' '014' '015' '016' '017' '018' '019' '020' ...
    '021' '022' '023' '024' '025' '026' '027' '028' '029' '030' ...
    '031' '032' '033' '034' '035' '101' '102' '103' '104' '105'} ;

% Define sea-ice edge concentration threshold
SIEdgeDefConc = 15 ;

% Load CESM mask for trimming
load('CESMmask.mat') % Ocean = 1; land = 0, per CESM grid

% Load XY points of LIA
load([XYPath 'XY_lia.mat'])
xstart = xstartlia ;
ystart = ystartlia ;

% Define time direction (forward or backward) and number of years of tracking
timearrow = -1 ; % backward tracking
nyears    = 10 ;

% Define the startmonth and all startyears for the tracking
startmonth = 3 ;
firstyear  = 1990 ; % 1990
lastyear   = 2100 ; % 2100

jstarts = YearMonth2JMonth_PDR(firstyear,startmonth):12:YearMonth2JMonth_PDR(lastyear,startmonth) ;

%% Eliminate starting points that are land in the CESM grid
maskstart = NaN(length(xstart),1) ;
for i = 1:length(xstart)
    maskstart(i) = CESMmask(xstart(i), ystart(i)) ;
end
b = find(maskstart == 0) ;
xstart(b) = [] ;
ystart(b) = [] ;
clear maskstart b

%% Looping through all the ensemble members
tic
for imem = 1:length(member)
    disp(['Member ' member{imem}])

    %% Looping through all the time slices
    for timeslice = 1:length(jstarts)
        jstart         = jstarts(timeslice) ; 
        jend           = jstart + timearrow*nyears*12 ; % +/- 1 * run time (months)
        [syear,smonth] = JMonth2YearMonth_PDR(jstart) ;
        [eyear,emonth] = JMonth2YearMonth_PDR(jend) ;
        
        disp(['Starting Year: ' int2str(syear)])

        % Declaring initial matrices
        % XLag & YLag: rows = individual grid cells in LIA; cols = months of run
        XLag            = NaN(numel(xstart),abs(jend-jstart)+1) ;
        YLag            = XLag ;
        XLag(:,1)       = xstart ;
        YLag(:,1)       = ystart ;
        ActiveFlag      = zeros(numel(xstart),2) ;
        ActiveFlag(:,1) = 1 ; % Initialize all as active

        %% Looping from jstart to jend +/- 1 (depending on time direction) 
        for jmonth = jstart:timearrow:jend-1*timearrow 
            n = abs(jmonth - jstart) + 1 ;
            
            % Find corresponding year and month for the interpolant
            if timearrow == 1
                [interpyear, interpmonth] = JMonth2YearMonth_PDR(jmonth) ;
            else
                [interpyear, interpmonth] = JMonth2YearMonth_PDR(jmonth-1) ;
            end
            
            % Load interpolant
            load([CESMLEInterpMPath member{imem} '/IceMotionInterpolant_CESMLE_' ...
                member{imem} '_' sprintf('%02d',interpmonth) '.' ...
                int2str(interpyear) '.mat']) ;
            
            % Determine time duration of the advection in s
            dt = eomday(interpyear,interpmonth)*86400 ;
            if eomday(interpyear,interpmonth) == 29
                dt = 28*86400 ; % No leap years in CESM
            end
            % The velocities are stored in cm/s
            fac = dt*.01/25000 ; % 25000 meters/EASEgrid cell and 100 cm/m
            clear interpyear interpmonth dt

            % Interpolate velocities from EASE grid onto exact location on that grid
            Ice = find(ActiveFlag(:,1)==1) ;
            if ~isempty(Ice) % Skip if no still-active ice floes
                uwork         = Fu(XLag(Ice,n),YLag(Ice,n)) ;
                vwork         = Fv(XLag(Ice,n),YLag(Ice,n)) ; 
                XLag(Ice,n+1) = XLag(Ice,n) + uwork*fac*timearrow ; 
                YLag(Ice,n+1) = YLag(Ice,n) + vwork*fac*timearrow ;
            end
            clear Fu Fv uwork vwork fac
            
            % Load SIC at jmonth after one time step
            [sicyear, sicmonth] = JMonth2YearMonth_PDR(jmonth+timearrow) ;
            load([CESMLEEASEGridMPath member{imem} '/sicEASE_CESMLE_' ...
                member{imem} '_' int2str(sicyear) '.mat'],'sicEASE') % only load sicEASE, not LAT & LON
            sic = sicEASE(:,:,sicmonth) ;
            clear sicEASE

            % Find the sea-ice concentration of the closest grid point and
            % check if over land in the CESM mask
            SIC  = NaN(length(ActiveFlag),1) ; % !!! Change for zeros instead?!?!
            mask = NaN(length(ActiveFlag),1) ;
            for i = 1:length(Ice)
                SIC(Ice(i))  = sic(round(XLag(Ice(i),n+1)),round(YLag(Ice(i),n+1))) ;
                mask(Ice(i)) = CESMmask(round(XLag(Ice(i),n+1)),round(YLag(Ice(i),n+1))) ;
            end

            % Find the buoys that were active and are still over ice
            ActiveFlag(:,2) = (SIC >= SIEdgeDefConc & ActiveFlag(:,1) == 1) ;
            clear SIC

            % Set the buoys that were active and have moved over land to -1
            ActiveFlag(mask == 0 & ActiveFlag(:,1) == 1,2) = -1 ;
            clear mask

            % Set the buoys that were active and have moved south of the 
            % Bering Strait to NaNs as well as their new position
            ActiveFlag(round(YLag(:,n+1)) < 79 & ActiveFlag(:,1)==1,2) = NaN ;
            XLag(isnan(ActiveFlag(:,2)) & ActiveFlag(:,1)==1,n+1) = NaN ;
            YLag(isnan(ActiveFlag(:,2)) & ActiveFlag(:,1)==1,n+1) = NaN ;
            
            % If the buoy has melted, set its new position to NaN so that
            % its previous position is the final one. This is because we
            % consider ice formation to occur at time t+1 when the grid
            % cell goes from open water at time t to ice at time t+1. Here,
            % because we are going backward in time, the final position
            % needs to match with when the ice formed.
            XLag(ActiveFlag(:,2)==0 & ActiveFlag(:,1)==1,n+1) = NaN ;
            YLag(ActiveFlag(:,2)==0 & ActiveFlag(:,1)==1,n+1) = NaN ;

            % Move parcels that are now over land to the last ocean grid
            % cell along their trajectory
            OnLand = find(ActiveFlag(:,2)==-1 & ActiveFlag(:,1)==1) ;
            if ~isempty(OnLand)
                for j = 1:length(OnLand)

                    % Declare temporary x and y variables and set to final
                    % position
                    xtemp = XLag(OnLand(j),n+1) ;
                    ytemp = YLag(OnLand(j),n+1) ;

                    % Calculate the angle made by a linear trajectory between
                    % initial and final positions in absolute value
                    angle = abs(atan((YLag(OnLand(j),n+1) - YLag(OnLand(j),n))/...
                        (XLag(OnLand(j),n+1) - XLag(OnLand(j),n)))) ;

                    % Determine the sign of delta x
                    if XLag(OnLand(j),n+1) > XLag(OnLand(j),n)
                        xsign = 1 ;
                    elseif XLag(OnLand(j),n+1) < XLag(OnLand(j),n)
                        xsign = -1 ;
                    else
                        xsign = 0 ;
                    end

                    % Determine the sign of delta y
                    if YLag(OnLand(j),n+1) > YLag(OnLand(j),n)
                        ysign = 1 ;
                    elseif YLag(OnLand(j),n+1) < YLag(OnLand(j),n)
                        ysign = -1 ;
                    else
                        ysign = 0 ;
                    end

                    % Calculate the value of CESM mask at final position 
                    % (should always be 0 since we haven't moved the parcel
                    % yet)
                    masktemp = CESMmask(round(xtemp),round(ytemp)) ;

                    % Loop until the parcel is no longer over a grid cell where
                    % mask is 0
                    while masktemp == 0

                        % Displace the parcel from 1 unit grid cell towards its
                        % initial position
                        xtemp = xtemp - xsign*cos(angle) ;
                        ytemp = ytemp - ysign*sin(angle) ;

                        % Make sure that the parcel hasn't been moved behind
                        % its initial position
                        nit = 0 ;
                        while (xtemp - XLag(OnLand(j),n))*xsign < 0 || ...
                                (ytemp - YLag(OnLand(j),n))*ysign < 0
                            xtemp = xtemp + 0.1*xsign*cos(angle) ;
                            ytemp = ytemp + 0.1*ysign*sin(angle) ;
                            nit = nit + 1 ;

                            % If the increment of 0.1 still brings the parcel
                            % behind its initial position, set to initial value
                            if nit == 10
                                xtemp = XLag(OnLand(j),n) ;
                                ytemp = YLag(OnLand(j),n) ;
                            end
                        end

                        % Calculate the value of CESM mask at the new position
                        masktemp = CESMmask(round(xtemp),round(ytemp)) ;

                        % If the parcel is now over a grid cell where CESM mask
                        % is 1, replace the original final position with the
                        % temporary position and change ActiveFlag to 1
                        if masktemp == 1
%                             disp('Safe!')
                            XLag(OnLand(j),n+1) = xtemp ;
                            YLag(OnLand(j),n+1) = ytemp ;
                            ActiveFlag(OnLand(j),2) = 1 ;
                        end
                    end
                end
            end
            clear OnLand xtemp ytemp angle masktemp

            % Update ActiveFlag col 1 with col 2  
            ActiveFlag(isnan(ActiveFlag(:,2)),2) = 0 ; % Set not active to 0
            ActiveFlag(:,1) = ActiveFlag(:,2) ;
            ActiveFlag(:,2) = 0 ; % Set flag of next position to 0

        end
        clear ActiveFlag

        %% Save data
        save(['Data_BackTrack_LIA_CESMLE_' member{imem} '_' int2str(syear) '-' ...
            sprintf('%02d',smonth) 'to' int2str(eyear) '-' sprintf('%02d',emonth) '.mat'],...
            'XLag','YLag')
        clear XLag YLag

        movefile(['Data_BackTrack_LIA_CESMLE_' member{imem} '_' int2str(syear) '-' ...
            sprintf('%02d',smonth) 'to' int2str(eyear) '-' sprintf('%02d',emonth) '.mat'],...
            CESMLELIAMPath)

    end
    toc
end
