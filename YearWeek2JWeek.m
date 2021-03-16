function [JWeek] = YearWeek2JWeek_PDR(Year,Week)
% YearWeek2JWeek a function to create the "Julian"
% week, with 1/1/1979 as the starting point.  
% Needed to work with weekly SSMI velocity data,
% where we need to create lagrangian tracks across
% year boundaries.  

refyear = 1979 ;
jyear = Year - refyear;
JWeek = 52 * jyear + Week ;
end

