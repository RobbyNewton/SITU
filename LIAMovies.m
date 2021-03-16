% LIAPlots
%
% Makes a movie of backtracks from the LIA.
% Selects 5 ensemble members to show, randomly.
% 
% BN, May 2019

FAITPaths_BN
TromsoPatches
load CESM_TracksInArrays.mat
load bathyEASE.mat
members = ['001';'002';'003';'004';'005';'006';'007';'008';'009';'010';...
    '011';'012';'013';'014';'015';'016';'017';'018';'019';'020';'021';...
    '022';'023';'024';'025';'026';'027';'028';'029';'030';'031';'032';...
    '033';'034';'035';'101';'102';'103';'104';'105'] ;
TimeSlices = {'2020_2010';'2030_2020';'2040_2030';'2050_2040'; ...
              '2060_2050';'2070_2060';'2080_2070';'2090_2080'} ;
jstarts = 2041 :120: 2881 ;
M = length(members) ;               % Number of ensemble members
T = length(TimeSlices) ;            % Number of time slices.
I = length(indLIA) ;                % Number of points in the LIA patch.
J = size(XLagsMnLIA2020_2010,2) ;   % Number of months in each time slice.

figure(2)
c = [.54 .5 .5] ;
figure(1), clf
mesh(bathyEASE,'FaceColor','Interp')
view(0,90),axis equal,hold on
fill(ab,bb,c,ac,bc,c,ae,be,c,ak,bk,c,al,bl,c,alia,blia,MyColors.IcyBlue, ...
    an,bn,c,ap,bp,c,az,bz,c) 
set(gca,'XLim',[70 265],'YLim',[90 280])
set(gca,'XTick',[],'YTick',[])
print -djpeg99 ArcticMap_TromsoPatches_EASEGrid.jpg

% Color schemes:
% % a = get(gcf,'colormap') ;
% % a(end,:) = [.4 .5 .3] ;
% % set(gcf,'colormap',a)
% % % or
colormap('bone')
a = get(gcf,'colormap') ;
a(end,:) = MyColors.EarthyGreen ; % .5 .35 0] ;     % Set the land to brown-ish.

set(gcf,'colormap',a)

%%
for t = 1:length(TimeSlices)
    figure(3)
    jstart      = jstarts(t) ;
    timeslice   = TimeSlices{t} ;
    YLags       = eval(['YLagsLIA' timeslice]) ;
    XLags       = eval(['XLagsLIA' timeslice]) ;
    perm        = randperm(40) ;
    mems        = perm(1:5) ;
    memstr      = [int2str(mems(1)) '_' int2str(mems(2)) '_' ...
        int2str(mems(3)) '_' int2str(mems(4)) '_' int2str(mems(5))] ;
    for i = 1:122
        [y,m] = JMonth2YearMonth_PDR(jstart - i) ;
        title(['Year: ' int2str(y) ';  Month: ' pad(int2str(m),2,'left','0')], ...
            'FontWeight','Bold','FontSize',18)
        plot(squeeze(YLags(:,i,mems(1))), squeeze(XLags(:,i,mems(1))), 'g.')
        plot(squeeze(YLags(:,i,mems(2))), squeeze(XLags(:,i,mems(2))), 'r.')
        plot(squeeze(YLags(:,i,mems(3))), squeeze(XLags(:,i,mems(3))), 'y.')
        plot(squeeze(YLags(:,i,mems(4))), squeeze(XLags(:,i,mems(4))), 'c.')
        plot(squeeze(YLags(:,i,mems(5))), squeeze(XLags(:,i,mems(5))), 'm.')
        p = get(gcf,'Position') ;
        LIAM(i) = getframe(gcf) ;
        unplot
        unplot
        unplot
        unplot
        unplot
    end
    v = VideoWriter(['LIABackTrack_' timeslice '_Mems_' memstr '.avi']) ;
    v.FrameRate = 5 ;
    open(v) ;
    for i = 1:length(LIAM)
        frame = LIAM(i) ;
        writeVideo(v,frame) ;
    end
    close(v)
end
