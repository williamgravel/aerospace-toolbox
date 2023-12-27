function [toTime] = convertTime(fromTime,fromScale,fromFormat,toScale,toFormat)
%CONVERTTIME Convert and format times between atomic time scales.
%   This function is capable of interpreting times expressed in various atomic time scales and
%   time formats and converting these to other time scales and time formats by using TAI time
%   as the common time scale (i.e. fromTime -> TAI -> toTime).
%
%   Inputs:
%
%      fromTime                - input time (1xN double/1x1 datetime)
%      fromScale               - input time scale (1x1 string)
%      fromFormat              - input time format (1x1 string)
%      toScale                 - output time scale (1x1 string)
%      toFormat                - output time format (1x1 string)
%
%   Outputs:
%
%      toTime                  - output time (1x1 double/datetime)
%
%   Notes:
%    - Allowable time scales comprise some of the common atomic time standards, including:
%       (a) 'tai'              - International Atomic Time
%       (b) 'utc'              - Coordinated Universal Time
%       (c) 'gps'              - GPS (USA)
%       (d) 'glonass'          - GLONASS (RUS)
%       (e) 'galileo'          - Galileo (EUR)
%       (f) 'beidou'           - BeiDou (CHN)
%       (g) 'qzss'             - QZSS/Michibiki (JPN)
%       (h) 'irnss'            - IRNSS/NavIC (IND)
%
%    - Allowable time formats include:
%       (i) 'tow'              - time of week vector = [TOW, (WEEK), (CYCLE)]
%       (j) 'dtobj'            - datetime object
%       (k) 'dtnum'            - datetime vector = [YEAR, MONTH, DAY, (HOUR, MINUTE, SECOND)]
%       (l) 'jd'               - Julian date
%       (m) 'mjd'              - modified Julian date
%
%    - GLONASS time is interpreted to be identical to UTC time, excluding the MSK timezone UTC+3
%      offset (this is because MATLAB does not provide a way for a local timezone to account for
%      leap seconds).
%
%    - All GNSS times except for GLONASS time can be express in TOW format. Also, datetime
%      objects for GPS, GALILEO, BEIDOU, QZSS, IRNSS, and TAI are formatted with the 'UTC'
%      timezone, while other time scales are represented with the 'UTCLeapSeconds' timezone.
%
%    - Input times expressed in TOW format can take the following forms:
%      (n) [t]                 - t = time since epoch
%      (o) [tow, week]         - tow = time of week; week = week number since epoch
%      (p) [tow, week, cycle]  - tow = time of week; week = week number; cycle = cycle number
%
%    - Input times expressed in DTNUM format can take the following forms:
%      (r) [Y, M, D]           - YMD = year, month, day
%      (q) [Y, M, D, H, MI, S] - YMD = year, month, day; HMIS = hour, minute, second
%
%   Author:     William Gravel
%   Created:    09/30/2022
%   Edited:     10/24/2022
%
%   See also DATETIME, JULIANDATE, SECONDS, MUSTBEVALIDTIMESTANDARD.

arguments
    fromTime (:,:) {mustBeA(fromTime,{'double','datetime'})}
    fromScale (1,1) string {mustBeValidTimeStandard}
    fromFormat (1,1) string {mustBeMember(fromFormat,{'tow','dtobj','dtnum','mjd','jd'})}
    toScale (1,1) string {mustBeValidTimeStandard}
    toFormat (1,1) string {mustBeMember(toFormat,{'tow','dtobj','dtnum','mjd','jd'})}
end

% Force input/output names to be lowercase
fromScale = lower(fromScale);
fromFormat = lower(fromFormat);
toScale = lower(toScale);
toFormat = lower(toFormat);

% Verify input/output time can be expressed in TOW if requested
if any(strcmp(fromScale,{'glonass','utc','tai'})) && strcmp(fromFormat,'tow')
    error('Unable to express %s time in "time of week" format.',upper(fromScale));
elseif any(strcmp(toScale,{'glonass','utc','tai'})) && strcmp(toFormat,'tow')
    error('Unable to express %s time in "time of week" format.',upper(toScale));
end

% Define cycle periods and epochs for TOW-based time systems
GNSS_CT = ["gps","galileo","beidou","qzss","irnss"];
GNSS_LS = "glonass";
GNSS_ALL = [GNSS_CT,GNSS_LS];

TS_CT = [GNSS_CT,"tai"];
TS_LS = [GNSS_LS,"utc"];
TS_ALL = [TS_CT,TS_LS];

W2S = 7*24*60*60;

C2W_GPS = 1024;
C2W_GAL = 4096;
C2W_BDT = 8192;
C2W_QZS = 1024;
C2W_IRN = 1024;

C2W = dictionary(GNSS_CT,[C2W_GPS,C2W_GAL,C2W_BDT,C2W_QZS,C2W_IRN]);

EPOCH_GPS = datetime(1980,1,6,0,0,0,'TimeZone','UTC');
EPOCH_GAL = datetime(1999,8,21,23,59,47,'TimeZone','UTC');
EPOCH_BDT = datetime(2006,1,1,0,0,0,'TimeZone','UTC');
EPOCH_QZS = datetime(1980,1,6,0,0,0,'TimeZone','UTC');
EPOCH_IRN = datetime(1999,8,21,23,59,47,'TimeZone','UTC');

EPOCH = dictionary(GNSS_CT,[EPOCH_GPS,EPOCH_GAL,EPOCH_BDT,EPOCH_QZS,EPOCH_IRN]);

TIMEZONE = dictionary(TS_ALL,["UTC","UTC","UTC","UTC","UTC","UTC","UTCLeapSeconds","UTCLeapSeconds"]);

%% Input Time Parser
switch fromFormat
    case 'tow'
        if size(fromTime,2) == 1
            tow = fromTime;
            week = zeros(size(tow));
            cycle = zeros(size(tow));
        elseif size(fromTime,2) == 2
            tow = fromTime(:,1);
            week = fromTime(:,2);
            cycle = zeros(size(tow));
        elseif size(fromTime,2) == 3
            tow = fromTime(:,1);
            week = fromTime(:,2);
            cycle = fromTime(:,3);
        else
            error('Incorrect number of fields.')
        end

        inputTime = EPOCH(fromScale) + seconds(tow + week*W2S + cycle*C2W(fromScale)*W2S);
    case 'dtnum'
        inputTime = datetime(fromTime,"TimeZone",TIMEZONE(fromScale));
    case 'dtobj'
        if isa(fromTime,"datetime") && ~strcmp(fromTime.TimeZone,TIMEZONE(fromScale))
            if strcmp(TIMEZONE(fromScale),"UTC")
                inputTime = datetime(fromTime,"TimeZone","UTC");
            elseif strcmp(TIMEZONE(fromScale),"UTCLeapSeconds")
                [y,m,d] = ymd(fromTime);
                [h,mi,s] = hms(fromTime);
                inputTime = datetime([y,m,d,h,mi,s],"TimeZone","UTCLeapSeconds");
            end
        end
    case 'jd'
        inputTime = datetime(fromTime,"ConvertFrom","juliandate","TimeZone",TIMEZONE(fromScale));
    case 'mjd'
        inputTime = datetime(fromTime,"ConvertFrom","modifiedjuliandate","TimeZone",TIMEZONE(fromScale));
end

%% Conversion Case Selection
if any(strcmp(fromScale,TS_LS)) && any(strcmp(toScale,TS_LS))                  % GLONASS/UTC -> GLONASS/UTC
    outputTime = inputTime;
elseif any(strcmp(fromScale,TS_LS))                                             % GLONASS/UTC -> GPS/GALILEO/BEIDOU/QZSS/IRNSS/TAI
    ls = leapseconds;
    outputTime = utc2tai(inputTime,ls) - TAI_OFFSET(fromScale);
elseif any(strcmp(toScale,TS_LS))                                               % GPS/GALILEO/BEIDOU/QZSS/IRNSS/TAI -> GLONASS/UTC
    ls = leapseconds;
    outputTime = tai2utc(inputTime,ls) + TAI_OFFSET(toScale);
else                                                                             % GPS/GALILEO/BEIDOU/QZSS/IRNSS/TAI -> GPS/GALILEO/BEIDOU/QZSS/IRNSS/TAI
    if ~strcmp(fromScale,"tai") && inputTime < EPOCH(fromScale)
        error('Input time before %s epoch is invalid.',upper(fromScale))
    end

    outputTime = inputTime + TAI_OFFSET(fromScale) - TAI_OFFSET(toScale);

    if ~strcmp(toScale,"tai") && outputTime < EPOCH(toScale)
        error('Output time before %s epoch is invalid.',upper(toScale))
    end
end

%% Output Time Formatter
switch toFormat
    case 'tow'
        tow = mod(floor(seconds(outputTime - EPOCH(toScale))),W2S);
        week = mod(floor(seconds(outputTime - EPOCH(toScale))/W2S),C2W(toScale));
        cycle = floor(seconds(outputTime - EPOCH(toScale))/W2S/C2W(toScale));
        
        toTime = [tow, week, cycle];
    case 'dtnum'
        [Y, M, D] = ymd(outputTime);
        [H, MI, S] = hms(outputTime);
        toTime = [Y, M, D, H, MI, S];
    case 'dtobj'
        toTime = outputTime;
    case 'jd'
        toTime = juliandate(outputTime,"juliandate");
    case 'mjd'
        toTime = juliandate(outputTime,"modifiedjuliandate");
end

%% Extra Information
% cycle = ((GNSST - tow)/W2S - week)/C2W;

% 1: datetime object
% 1: [tow, week, cycle]
% 1: [tow, week] (assume current cycle)
% 6: y, m, d, h, mi, s
% 3: y, m, d
% 1: mjd

% GPS %
% Epoch:  1980-01-06 00:00:00 (UTC)
% 1st RO: 1999-08-22 00:00:00 (GPST), 1999-08-21 23:59:47 (UTC)
% 2nd RO: 2019-04-07 00:00:00 (GPST), 2019-04-06 23:59:42 (UTC)
% 3rd RO: 2038-11-21 00:00:00 (GPST0, TBD (UTC)

% Galileo %
% Epoch:  1999-08-22 00:00:00 (GPST), 1999-08-21 23:59:47 (UTC)
% 1st RO: 2078-02-20 00:00:00 (GPST), TBD (UTC)

% GLONASS %
% UTC + 3h (except in RINEX, no offset, i.e. original UTC)

% Beidou %
% Epoch:  2006-01-01 00:00:14 (GPST), 2006-01-01 00:00:00 (UTC)
% 1st RO: 2161-02-01 00:00:14 (GPST), TBD (UTC)

% QZSS %
% Epoch:  1999-08-22 00:00:00 (GPST), 

% IRNSS %
% Epoch:  1999-08-22 00:00:00 (GPST), 1999-08-21 23:59:47 (UTC)
% 1st RO: 2019-04-07 00:00:00 (GPST), 2019-04-06 23:59:42 (UTC)
% 2nd RO: 2038-11-21 00:00:00 (GPST), TBD (UTC)

end
