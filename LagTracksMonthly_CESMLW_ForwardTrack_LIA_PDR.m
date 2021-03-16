%  LagTrackMonthly_CESMLW_ForwardTrack_LIA_PDR
%
%  A scrip to:
%   * Load the coordinates of the newly formed ice found using
%   EEZ_CESMLW_NewIce_Monthly_PDR.m
%   * Track each ice parcel from the moment it was formed to a maximum
%   number of years
%   * Save all positions of tracks that enter the LIA
%
%  This script is a modified version of LagTracksEEZMonthly_CESMLW_PDR.m
%  that was used to track ice parcels forward in time for the EEZFuture
%  project (Earth's Future (2020) paper).
%
%  Patricia, June 2020
%  Last update: July 2020

%% Common Variables
% % When working on Mac
% addpath('/Users/pade7652/Google Drive/CUBoulder/Projects/LIA/Patricia/Scripts/')
% FAITPaths_Mac_LIA_PDR

% When working on Cheyenne
FAITPaths_Cheyenne_LIA

% Read EASE Grid coordinates
load('north_x_y_lat_lon')
x       = north_x_y_lat_lon(:,1) + 1 ;
y       = north_x_y_lat_lon(:,2) + 1 ;
XX      = reshape(x,361,361) ;
YY      = reshape(y,361,361) ;
clear north_x_y_lat_lon

% Ensemble members
member = {'001' '002' '003' '004' '005' '006' '007' '008' '009' '010' '011'} ;

% Define sea-ice edge concentration threshold
SIEdgeDefConc = 15 ;

% Maximum number of years of tracking
nyears = 20 ;

% Load CESM mask for trimming
load('CESMmask.mat') % Ocean = 1; land = 0, per CESM grid

% Load XY points of LIA and eliminate those that are over land
load([XYPath 'XY_lia.mat'])
maskstart = NaN(length(xstartlia),1) ;
for i = 1:length(xstartlia)
    maskstart(i) = CESMmask(xstartlia(i),ystartlia(i)) ;
end
b = find(maskstart == 0) ;
xstartlia(b) = [] ;
ystartlia(b) = [] ;
clear maskstart b

%% Looping through all the ensemble members
tic
for imem = 1:length(member)
    disp(['Member ' member{imem}])
 
    %% Load initial points for EEZ tracking - this loads xstart, ystart and jmonths
    
    % When working on Cheyenne
    load([CESMLWEEZMPath 'Formation/Matrix_NewIce_CESMLW_' member{imem} '.mat'])
    
%     % When working on Mac
%     load(['/Users/pade7652/Google Drive/CUBoulder/Projects/EEZFuture/Data/' ...
%         'CESMLW/monthly/Formation/Matrix_NewIce_CESMLW_' member{imem} '.mat'])

    %% Declaring initial matrices
    XLag       = NaN(numel(xstart),nyears*12+1) ; % maximum number of years of tracking + initial position
    YLag       = XLag ;
    time       = XLag ;
    XLag(:,1)  = xstart ;
    YLag(:,1)  = ystart ;
    
    XTemp      = NaN(numel(xstart),2) ;
    YTemp      = XTemp ;
    XTemp(:,1) = xstart ;
    YTemp(:,1) = ystart ;
    ActiveFlag = zeros(numel(xstart),2) ;
    column     = ones(numel(xstart),1) ;

    %% Looping through the time period
    jmonthsu = unique(jmonths) ;
    for jmonth = jmonthsu(1):1:jmonthsu(end) % Loop through all months
        
        % Find corresponding year and month
        [interpyear,interpmonth] = JMonth2YearMonth_PDR(jmonth) ;
        display(['Year: ' int2str(interpyear) ', Month: ' int2str(interpmonth)])
        
        % Load interpolant
        if interpyear <= 2005
            load([CESMLEInterpMPath member{imem} '/IceMotionInterpolant_CESMLE_' ...
                member{imem} '_' sprintf('%02d',interpmonth) '.' ...
                int2str(interpyear) '.mat']) ;
        else
            load([CESMLWInterpMPath member{imem} '/IceMotionInterpolant_CESMLW_' ...
                member{imem} '_' sprintf('%02d',interpmonth) '.' ...
                int2str(interpyear) '.mat']) ;
        end
        
        % Determine time duration of the advection in s
        dt = eomday(interpyear,interpmonth)*86400 ;
        if eomday(interpyear,interpmonth) == 29
            dt = 28*86400 ; % No leap years in CESM
        end
        % The velocities are stored in cm/s
        fac = dt*.01/25000 ; % 25000 meters/EASEgrid cell and 100 cm/m ;
        clear interpyear interpmonth dt

        % Activate the newly formed ice parcels and record the jmonth at
        % the time of formation
        ActiveFlag(jmonths == jmonth,1) = 1 ;
        time(jmonths == jmonth) = jmonth ;

        % Find the index of all active ice floes
        Ice = find(ActiveFlag(:,1) == 1) ;
        if ~isempty(Ice) % Skip if no active ice floes
            
            % Interpolate velocities from EASE grid onto exact location on that grid
            uwork        = Fu(XTemp(Ice,1),YTemp(Ice,1)) ;
            vwork        = Fv(XTemp(Ice,1),YTemp(Ice,1)) ;
            
            % Move active parcels to their new position after one time step
            XTemp(Ice,2) = XTemp(Ice,1) + uwork*fac ; 
            YTemp(Ice,2) = YTemp(Ice,1) + vwork*fac ;
            clear Fu Fv uwork vwork fac
            
            % Load SIC
            [sicyear,sicmonth] = JMonth2YearMonth_PDR(jmonth+1) ;
            if sicyear <= 2005
                load([CESMLEEASEGridMPath member{imem} '/sicEASE_CESMLE_' ...
                    member{imem} '_' int2str(sicyear) '.mat'])
            else
                load([CESMLWEASEGridMPath member{imem} '/sicEASE_CESMLW_' ...
                    member{imem} '_' int2str(sicyear) '.mat'])
            end
            sic = sicEASE(:,:,sicmonth) ;
            clear sicEASE

            % Find the sea-ice concentration of the closest grid point and
            % check if over NaN in the SIC or SIV field
            SIC = NaN(length(ActiveFlag),1) ;
            mask = NaN(length(ActiveFlag),1) ;
            for i = 1:length(Ice)
                SIC(Ice(i)) = sic(round(XTemp(Ice(i),2)),round(YTemp(Ice(i),2))) ;
                mask(Ice(i)) = CESMmask(round(XTemp(Ice(i),2)),round(YTemp(Ice(i),2))) ;
            end
            
            % Find the buoys that were active and are still over ice
            ActiveFlag(:,2) = (SIC >= SIEdgeDefConc & ActiveFlag(:,1) == 1) ;
            clear SIC

            % Set the buoys that were active and have moved over land to -1
            ActiveFlag(mask == 0 & ActiveFlag(:,1) == 1,2) = -1 ;
            clear mask

            % Set the buoys that were active and have moved south of the 
            % Bering Strait to NaNs as well as their new position
            ActiveFlag(round(YTemp(:,2)) < 79 & ActiveFlag(:,1) == 1,2) = NaN ;
            XTemp(isnan(ActiveFlag(:,2)) & ActiveFlag(:,1) == 1,2) = NaN ;
            YTemp(isnan(ActiveFlag(:,2)) & ActiveFlag(:,1) == 1,2) = NaN ;

            % Move parcels that are now over land to the last ocean grid
            % cell along their trajectory
            OnLand = find(ActiveFlag(:,2) == -1 & ActiveFlag(:,1) == 1) ;
            if ~isempty(OnLand)

                for j = 1:length(OnLand)

                    % Declare temporary x and y variables and set to final
                    % position
                    xtemp = XTemp(OnLand(j),2) ;
                    ytemp = YTemp(OnLand(j),2) ;

                    % Calculate the angle made by a linear trajectory between
                    % initial and final positions in absolute value
                    angle = abs(atan((YTemp(OnLand(j),2) - YTemp(OnLand(j),1))/...
                        (XTemp(OnLand(j),2) - XTemp(OnLand(j),1)))) ;

                    % Determine the sign of delta x
                    if XTemp(OnLand(j),2) > XTemp(OnLand(j),1) ; xsign = 1 ;
                    elseif XTemp(OnLand(j),2) < XTemp(OnLand(j),1) ; xsign = -1 ;
                    else ; xsign = 0 ;
                    end

                    % Determine the sign of delta y
                    if YTemp(OnLand(j),2) > YTemp(OnLand(j),1) ; ysign = 1 ;
                    elseif YTemp(OnLand(j),2) < YTemp(OnLand(j),1) ; ysign = -1 ;
                    else ; ysign = 0 ;
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
                        n = 0 ;
                        while (xtemp - XTemp(OnLand(j),1))*xsign < 0 || ...
                                (ytemp - YTemp(OnLand(j),1))*ysign < 0
                            xtemp = xtemp + 0.1*xsign*cos(angle) ;
                            ytemp = ytemp + 0.1*ysign*sin(angle) ;
                            n = n + 1 ;

                            % If the increment of 0.1 still bring the parcel
                            % behind its initial position, set to initial value
                            if n == 10
                                xtemp = XTemp(OnLand(j),1) ;
                                ytemp = YTemp(OnLand(j),1) ;
                            end
                        end

                        % Calculate the value of CESM mask at the new position
                        masktemp = CESMmask(round(xtemp),round(ytemp)) ;

                        % If the parcel is now over a grid cell where CESM mask
                        % is 1, replace the original final position with the
                        % temporary position and change ActiveFlag to 1
                        if masktemp == 1
    %                         disp('Safe!')
                            XTemp(OnLand(j),2) = xtemp ;
                            YTemp(OnLand(j),2) = ytemp ;
                            ActiveFlag(OnLand(j),2) = 1 ;
                        end
                    end
                end
            end
            clear OnLand xtemp ytemp angle masktemp
            
            % Find active buoy tracks that haven't reach the maximum number
            % of years and the ones that have
            LessYrs = find(ActiveFlag(:,1) == 1 & column < nyears*12+1 & ~isnan(ActiveFlag(:,2))) ;
            MaxYrs  = find(ActiveFlag(:,1) == 1 & column >= nyears*12+1) ;
            if ~isempty(MaxYrs)
                disp([int2str(length(MaxYrs)) ' track(s) have reached the maximum number of years.'])
            end

            % For LessYrs buoys:
            % Increase the column index by 1
            column(LessYrs) = column(LessYrs) + 1 ;
            
            % Record new position and time
            for k = 1:length(LessYrs)
                XLag(LessYrs(k),column(LessYrs(k))) = XTemp(LessYrs(k),2) ;
                YLag(LessYrs(k),column(LessYrs(k))) = YTemp(LessYrs(k),2) ;
                time(LessYrs(k),column(LessYrs(k))) = jmonth + 1 ;
            end
            
            % For MaxYrs buoys:
            % Set ActiveFlag(:,2) to 0
            ActiveFlag(MaxYrs,2) = 0 ;
            clear LessYrs MaxYrs
            
            % Copy column 2 into column 1 for next iteration
            XTemp(Ice,1) = XTemp(Ice,2) ;
            YTemp(Ice,1) = YTemp(Ice,2) ;
            XTemp(Ice,2) = NaN ; % Set next positions to NaN
            YTemp(Ice,2) = NaN ; 
            ActiveFlag(isnan(ActiveFlag(:,2)),2) = 0 ; % Set NaN values to 0
            ActiveFlag(:,1) = ActiveFlag(:,2) ;
            ActiveFlag(:,2) = 0 ; % Set flag of next position to 0

        end

    end % done looping through all months
    clear XTemp YTemp ActiveFlag column

    %% Removing that tracks that never entered the LIA

    insideLIA = zeros(size(XLag)) ;
    % Loop through the grid cells of the LIA
    for m = 1:length(xstartlia)
        insideLIA(round(XLag) == xstartlia(m) & round(YLag) == ystartlia(m)) = 1 ;
    end

    % Subset XLag, YLag and time with only the tracks that entered the LIA,
    % i.e. all tracks that have a value > 0 after summing over all columns
    XLag      = XLag(sum(insideLIA,2) > 0,:) ;
    YLag      = YLag(sum(insideLIA,2) > 0,:) ;
    time      = time(sum(insideLIA,2) > 0,:) ;
    insideLIA = insideLIA(sum(insideLIA,2) > 0,:) ;

    %% Save data
    save(['Data_ForwardTrack_' int2str(nyears) 'yrs_LIA_CESMLW_' member{imem} ...
        '.mat'],'XLag','YLag','time','insideLIA')
    clear XLag YLag time insideLIA

    movefile(['Data_ForwardTrack_' int2str(nyears) 'yrs_LIA_CESMLW_' member{imem} ...
        '.mat'],CESMLWLIAMPath)
    toc
    
end % done looping through all ensemble members
