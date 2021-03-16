function [JMonth] = YearMonth2JMonth_PDR( Year, Month )

% YearMonth2JMonth is a function to create the "Julian"
% month, with 1/1/1979 as the starting point.  
%
% Needed to work with monthly CCSM velocity data,
% where we need to create lagrangian tracks across
% year boundaries.  

refyear = 1850 ;
jyear = Year - refyear;
JMonth = 12 * jyear + Month ;
end

