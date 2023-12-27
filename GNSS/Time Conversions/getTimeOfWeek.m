function [TOW] = getTimeOfWeek(date)
%GETTIMEOFWEEK Determines time of week from given date.
%   This function calculates the time of week of a given a date and time for GNSS time scales. 
%   The time of week is calculated by computing the date's weekday and simply counting the
%   number of elapsed seconds within the week. The date's weekday is determined using the
%   Doomsday algorithm, which takes advantage of certain dates of the year sharing a common
%   weekday.
%
%   This program is especially useful since it does not require any datetime object conversions
%   (which are computationally expensive) or any information on the satellite time system in use
%   (since all GNSS time scales have the same TOW definitions).
%
%   Inputs:
%
%      date                    - full date and time [Y, M, D, H, MI, S] (1x6 double)
%
%   Outputs:
%
%      TOW                     - time of week (1x1 double)
%
%   Author:     William Gravel
%   Created:    10/24/2022
%   Edited:     10/25/2022
%
%   References:
%
%      Conway, J. H. (1973). Tomorrow is the day after Doomsday. Eureka, 36, 28–31.
%       https://archim.org.uk/eureka/archive/Eureka-36.pdf
%
%   See also CONVERTTIME.

arguments
    date (1,6) double
end

% Define anchor days and doomsdays
ANCHORS = [5, 3, 2, 0];
DOOMSDAYS = [4, 29, 7, 4, 9, 6, 11, 8, 5, 10, 7, 12;
    3, 28, 7, 4, 9, 6, 11, 8, 5, 10, 7, 12];

% Extract datetime components
Y = date(1);
M = date(2);
D = date(3);
H = date(4);
MI = date(5);
S = date(6);

% Determine year's doomsday
doomsday_year = floor(mod(Y,100)/12) + mod(mod(Y,100),12) + floor(mod(mod(Y,100),12)/4);

% Determine century's doomsday
doomsday_cent = ANCHORS(floor(mod(Y+200,400)/100) + 1);

% Determine doomsday
doomsday = mod(doomsday_year + doomsday_cent,7);

% Determine whether given year is leap year
is_leap_year = (mod(Y,4) == 0) && ((mod(Y,100) ~= 0) || (mod(Y,400) == 0));

% Calculate offset from doomsday
doomsday_offset = rem(DOOMSDAYS(~is_leap_year + 1,M) - D,7);

% Determine weekday
weekday = mod(doomsday - doomsday_offset,7);

% Convert date into time of week
TOW = weekday*24*60*60 + H*60*60 + MI*60 + S;

end
