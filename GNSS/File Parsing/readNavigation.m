function [ephem,iono,time,leap,file] = readNavigation(filename,opts)
%READNAVIGATION Read and parse RINEX navigation files.
%   WRITE DESCRIPTION HERE
%
%   Inputs:
%
%      filename                - navigation file name (1x1 string)
%      opts.OutputType         - output type (1x1 string)
%      opts.CleanOutput        - clean output (1x1 logical)
%
%   Outputs:
%
%      ephem                  - ephemeris data (MxN double/1x1 struct/MxN table)
%      iono                   - ionospheric correction constants (1x1 struct)
%      time                   - time system correction constants (1x1 struct)
%      file                   - file information (1x1 struct)
%
%   RINEX Notes:
%    - GPS/Galileo/BeiDou/QZSS/IRNSS entries follow similar formats based on Keplerian orbital
%      elements, while GLONASS/SBAS entries follow similar formats based on Cartesian state
%      elements.
%
%    - Time of clock (TOC) entry data reported by RINEX only includes last 2 digits of year.
%      Program assumes 2-digit years greater than current 2-digit year are in previous century,
%      while 2-digit years equal or less than current 2-digit year are in current century.
%    * Example: If current year is 2022, 2-digit years 23-99 are assumed to be 1923-1999, while
%      2-digit years 0-22 are assumed to be 2000-2022.
%
%    - Time tags for GLONASS reported by RINEX are not offset from UTC by 3 hours like the real
%      GLONASS time scale (i.e. GLO = UTC). Time tags for SBAS reported by RINEX are expressed
%      in GPS time frame.
%
%    - RINEX version 3.05 introduced a fifth data line for GLONASS entries containing status
%      flags, health flags, URA, and group delay difference. GLONASS entries in prior RINEX
%      versions only conclude 4 data lines.
%
%    - Week numbers reported by RINEX are:
%       (1) continuous without roll-overs (no modulo operation is applied),
%       (2) expressed in each entries' respective GNSS time scales,
%       (3) zero-aligned with the GPS epoch (except for BeiDou week numbers which keep their
%           BDT epoch reference)
%    * Example: Week number for Galileo on Oct. 21, 2022 is actually equal to 1208 but reported
%      as 2232 = 1208 + 1024 (number of weeks difference between GPS and Galileo epochs).
%
%   Program Notes:
%    - This program is capable of reading the following RINEX navigation files:
%      +--------------------+----------+----------+
%      | GNSS Constellation | RINEX V2 | RINEX V3 |
%      +--------------------+----------+----------+
%      | GPS                | *.yyn    | *_GN.rnx |
%      | GLONASS            | *.yyg    | *_RN.rnx |
%      | Galileo            |          | *_EN.rnx |
%      | BeiDou             |          | *_CN.rnx |
%      | QZSS               |          | *_JN.rnx |
%      | IRNSS              |          | *_IN.rnx |
%      | SBAS               | *.yyh    | *_SN.rnx |
%      | Mixed              |          | *_MN.rnx |
%      +--------------------+----------+----------+
%
%    - This program uses standardized notation for common parameters in order to maintain
%      consistency. The notation used by this program and each GNSS ICD is shown below:
%      +--------------------------+---------+------+---------+--------+------+--------+
%      |        Parameter         | Program | GPS  | Galileo | BeiDou | QZSS | IRNSS  |
%      +--------------------------+---------+------+---------+--------+------+--------+
%      | Issue of data, ephemeris | IODE    | IODE | IODNav¹ | AODE   | IODE | IODEC² |
%      | Issue of data, clock     | IODC    | IODC | IODNav¹ | AODC   | IODC | IODEC² |
%      | User range accuracy      | URA     | URA  | SISA³   | URA    | URA  | URA    |
%      | Total group delay        | T_GD    | T_GD | BGD     | T_GD   | T_GD | T_GD   |
%      | Clock corrections        | a_fn    | a_fn | a_fn    | a_n    | a_fn | a_fn   |
%      +--------------------------+---------+------+---------+--------+------+--------+
%      ¹ IODNav reported by Galileo is the IOD for the ephemeris corrections, clock corrections,
%        and signal-in-space accuracy.
%      ² IODEC reported by IRNSS is the IOD for the ephemeris and clock corrections.
%      ³ SISA reported by Galileo is the signal-in-space accuracy for a given SV and should not
%        be confouded for a user range accuracy (URA) as reported by the other constellations.
%
%   Author:     William Gravel
%   Created:    09/30/2022
%   Edited:     11/05/2022
%
%   See also READOBSERVATION, READALMANAC.

arguments
    filename (1,1) string
    opts.OutputType (1,1) string {mustBeMember(opts.OutputType,["array","struct","table"])} = "struct"
    opts.CleanOutput (1,1) logical = false
end

%% Setup
% Define number of variables and number of GNSS
N_VARS_PER_GNSS = [31,20,30,29,31,28,18];
N_GNSS = 7;

% Define GNSS names and helper functions
GNSS_FULL = ["GPS","GLONASS","Galileo","BeiDou","QZSS","IRNSS","SBAS"];
GNSS_PART = ["GPS","GLO","GAL","BDS","QZS","IRN","SBAS"];
GNSS_INIT = ["G","R","E","C","J","I","S"];
gnss_p2f = @(part) GNSS_FULL(contains(GNSS_PART,part));
gnss_i2f = @(init) GNSS_FULL(contains(GNSS_INIT,init));
gnss_i2n = @(init) find(contains(GNSS_INIT,init));

% Define time system names and helper functions
TS_FULL = ["GPS","GLONASS","Galileo","BeiDou","QZSS","IRNSS","SBAS","UTC"];
TS_PART = ["GP","GL","GA","BD","QZ","IR","SB","UT"];
ts_p2f = @(part) TS_FULL(contains(TS_PART,part));

% Define UTC identifiers
UTC_FULL = ["UTC","UTC(NIST)","UTC(USNO)","UTC(SU)","UTC(BIPM)","UTC(Europe Lab)","UTC(CRL)","UTC(NTSC)","UTC"];
UTC_PART = ["","NIST","USNO","SU","BIPM","Europe Lab","CRL","NTSC",""];
UTC_IDS = [0,1,2,3,4,5,6,7,8];
utc_n2f = @(n) UTC_FULL(find(n >= UTC_IDS,1,'last'));
utc_n2p = @(n) UTC_PART(find(n >= UTC_IDS,1,'last'));

% Define variable names and units
GNSS_INIT_A = ["G","E","C","J","I"];
GNSS_INIT_B = ["R","S"];

VARIABLE_NAMES_A = ["PRN","health","M_0","e","sqrt_A","Omega_0","i_0","omega","Delta_n","Omega_dot","i_0_dot","Crc","Crs","Cic","Cis","Cuc","Cus","a_f0","a_f1","a_f2","IODE","IODC","TOE","TOC","TOT","WN","URA"];
VARIABLE_NAMES_B = ["PRN","health","sat_x","sat_y","sat_z","sat_u","sat_v","sat_w","sat_u_dot","sat_v_dot","sat_w_dot"];
VARIABLE_NAMES = struct(...
    "GPS",[VARIABLE_NAMES_A,"T_GD","CFI","L2_codes","L2_P_flag"],...
    "GLONASS",[VARIABLE_NAMES_B,"tau_n","gamma_n","t_k","AOII","freq_num","T_GD","URA","status_flags","health_flags"],...
    "Galileo",[VARIABLE_NAMES_A,"T_GD_a","T_GD_b","data_sources"],...
    "BeiDou",[VARIABLE_NAMES_A,"T_GD_1","T_GD_2"],...
    "QZSS",[VARIABLE_NAMES_A,"T_GD","CFI","L2_codes","L2_P_flag"],...
    "IRNSS",[VARIABLE_NAMES_A,"T_GD"],...
    "SBAS",[VARIABLE_NAMES_B,"a_f0","a_f1","IODE","IODC","TOT","WN","URA"]);

VARIABLE_UNITS_A = ["","","rad","","m^(1/2)","rad","rad","rad","rad/s","rad/s","rad/s","m","m","rad","rad","rad","rad","s","s/s","s/s^2","","","s","s","s","wk","m"];
VARIABLE_UNITS_B = ["","","m","m","m","m/s","m/s","m/s","m/s^2","m/s^2","m/s^2"];
VARIABLE_UNITS = struct(...
    "GPS",[VARIABLE_UNITS_A,"s","h","",""],...
    "GLONASS",[VARIABLE_UNITS_B,"s","s/s","s","day","","s","m","",""],...
    "Galileo",[VARIABLE_UNITS_A,"s","s",""],...
    "BeiDou",[VARIABLE_UNITS_A,"s","s"],...
    "QZSS",[VARIABLE_UNITS_A,"s","","",""],...
    "IRNSS",[VARIABLE_UNITS_A,"s"],...
    "SBAS",[VARIABLE_UNITS_B,"s","s/s","","","s","wk","m"]);

% Determine current year for full-year estimations in time tags (necessary for RINEX 2.XX)
TODAY_Y = year(datetime('today'));

%% Header Section
% Open file
lines = readlines(filename,'EmptyLineRule','skip');

% Read and parse file information depending on version number format
if lines{1}(7) == '.' % float -> version = X.XX
    finfo = textscan(lines{1}(1:60),'%9.2f %*11c %1c %*19c %1c %*19c','Whitespace','');
else % integer -> version = X
    finfo = textscan(lines{1}(1:60),'%6f %*14c %1c %*19c %1c %*19c','Whitespace','');
end

% Determine file version
fversion = finfo{1};
fversion_M = floor(fversion);

% Extract file constellation depending on RINEX version
if fversion_M == 2
    ftype = 'N';
    switch finfo{2}
        case 'N'
            fgnss_i = 'G';
        case 'G'
            fgnss_i = 'R';
        case 'H'
            fgnss_i = 'S';
    end
elseif fversion_M == 3
    ftype = finfo{2};
    fgnss_i = finfo{3};
else
    exception = MException('readNavigation:invalidFileVersion','Invalid file version number. Must be RINEX version 2 or 3.');
    throw(exception)
end
fgnss_n = gnss_i2n(fgnss_i);
fgnss_f = gnss_i2f(fgnss_i);

% Validate file type (must be navigation file)
if ftype ~= 'N'
    exception = MException('readNavigation:invalidFileType','Invalid file type. Must be RINEX navigation file.');
    throw(exception)
end

% Define number of lines per entry for each GNSS
if fversion < 3.05
    N_LINES_PER_ENTRY = [8,4,8,8,8,8,4];
else
    N_LINES_PER_ENTRY = [8,5,8,8,8,8,4];
end

% Assemble file information struct
file = struct('version',fversion,'type',ftype,'gnss',fgnss_i);

% Initialize ionospheric and time correction structs
iono = struct;
time = struct;
leap = struct;

% Read and parse remaining header information
N_headers = 1;
cline = '';
while ~contains(cline,'END OF HEADER')
    cline = lines(N_headers);
    cvalue = cline(1:60);
    clabel = cline(61:end);
    N_headers = N_headers + 1;

    if fversion_M == 2
        if contains(clabel,'ION ALPHA')
            pvalue = textscan(cvalue,'%12.4f',4,'Whitespace','');

            iono.GPS.alpha = [pvalue{:}]';
        elseif contains(clabel,'ION BETA')
            pvalue = textscan(cvalue,'%12.4f',4,'Whitespace','');

            iono.GPS.beta = [pvalue{:}]';
        elseif contains(clabel,'DELTA-UTC: A0,A1,T,W')
            pvalue = textscan(cvalue,'%*3c %19.12f %19.12f %9f %9f');

            time.GPS.ts_ref = 'UTC';
            time.GPS.A_0 = pvalue{1};
            time.GPS.A_1 = pvalue{2};
            time.GPS.TOW_ref = pvalue{3};
            time.GPS.WN_ref = pvalue{4};
        elseif contains(clabel,'CORR TO SYSTEM TIME')
            if fgnss_i == 'R'
                pvalue = textscan(cvalue,'%6f %6f %6f %*3c %19.12f');
        
                time.GLONASS.ts_ref = 'UTC(SU)';
                time.GLONASS.t_ref = datetime([pvalue{1:3}]);
                time.GLONASS.tau_C = -pvalue{4};
            elseif fgnss_i == 'S' && fversion < 2.11
                pvalue = textscan(cvalue,'%6f %6f %6f %*3c %19.12f');
        
                time.SBAS.ts_ref = 'UTC';
                time.SBAS.t_ref = datetime([pvalue{1:3}]);
                time.SBAS.W_0 = pvalue{4};
            end
        elseif contains(clabel,'D-UTC A0,A1,T,W,S,U')
            pvalue = textscan(cvalue,'%19.12f %19.12f %7f %5f %*1c %5c %*1c %2f %*2c');

            time.SBAS.ts_ref = utc_n2f(pvalue{6});
            time.SBAS.A_0 = pvalue{1};
            time.SBAS.A_1 = pvalue{2};
            time.SBAS.TOW_ref = pvalue{3};
            time.SBAS.WN_ref = pvalue{4};
            time.SBAS.sys = pvalue{5};
        elseif contains(clabel,'LEAP SECONDS')
            pvalue = textscan(cvalue,'%6f','Whitespace','');
            
            leap.count = pvalue{1};
        end
    elseif fversion_M == 3
        if contains(clabel,'IONOSPHERIC CORR')
            pvalue = textscan(cvalue,'%3c %1c %*1c %12.4f %12.4f %12.4f %12.4f %*1c %1c %*1c %2f','Whitespace','');
    
            egnss_f = gnss_p2f(pvalue{1});
            switch pvalue{2}
                case 'A'
                    iono.(egnss_f).alpha = [pvalue{3:6}];
                    iono.(egnss_f).alpha_TOT = [pvalue{7}];
                    iono.(egnss_f).alpha_sv = [pvalue{8}];
                case 'B'
                    iono.(egnss_f).beta = [pvalue{3:6}];
                    iono.(egnss_f).beta_TOT = [pvalue{7}];
                    iono.(egnss_f).beta_sv = [pvalue{8}];
                otherwise
                    iono.(egnss_f).ai = [pvalue{3:5}];
                    iono.(egnss_f).ai_TOT = [pvalue{7}];
                    iono.(egnss_f).ai_sv = [pvalue{8}];
                    
            end
        elseif contains(clabel,'TIME SYSTEM CORR')
            pvalue = textscan(cvalue,'%2c %2c %*1c %17.10f %16.9f %*1c %6f %*1c %4f %*1c %2c %1c %2f %*1c %2f %*1c','Whitespace','');
    
            ts_obs = ts_p2f(pvalue{1});
            ts_ref = ts_p2f(pvalue{2});
    
            time.(ts_obs).ts_ref = ts_ref;
            if fgnss_i == 'R'
                time.(ts_obs).tau_C = -pvalue{3};
            else
                time.(ts_obs).A_0 = pvalue{3};
                time.(ts_obs).A_1 = pvalue{4};
                time.(ts_obs).TOW_ref = pvalue{5};
                time.(ts_obs).WN_ref = pvalue{6};
            end
            time.(ts_obs).sid = pvalue{7};
            time.(ts_obs).sys = pvalue{8};
            time.(ts_obs).prn = pvalue{9};
            time.(ts_obs).utc = pvalue{10};
        elseif contains(clabel,'LEAP SECONDS')
            pvalue = textscan(cvalue,'%6f %6f %6f %6f %3c','Whitespace','');
    
            leap.count = pvalue{1};
            leap.future = pvalue{2};
            leap.future_WN = pvalue{3};
            leap.future_DN = pvalue{4};
            leap.ts = pvalue{5};
        end
    end
end

%% Setup
% Pre-allocate memory for each GNSS ephemerides by tallying up respective number of entries and
% initializing ephemeris arrays.
fid_count = fopen(filename,'r');
fdata = fread(fid_count);
fclose(fid_count);

N_lines = sum(fdata == 10);
N_entries_per_gnss = zeros(1,N_GNSS);

if fgnss_i == 'M'
    feol_ind = find(fdata == 10);
    egnss_all = fdata(feol_ind(N_headers:end-1) + 1);

    for egnss_n = 1:N_GNSS
        N_entries_per_gnss(egnss_n) = sum(egnss_all == double(GNSS_INIT(egnss_n)));
    end
else
    N_entries_per_gnss(fgnss_n) = (N_lines - N_headers)/N_LINES_PER_ENTRY(fgnss_n);
end

ephem = struct;
for egnss_n = 1:N_GNSS
    if N_entries_per_gnss(egnss_n) ~= 0
        ephem.(GNSS_FULL(egnss_n)) = zeros(N_entries_per_gnss(egnss_n),N_VARS_PER_GNSS(egnss_n));
    end
end

%% Data Record
% Initialize entry count for each GNSS
n_line = ones(1,N_GNSS);

while true
    if feof(fid)
        break
    end

    % Interpret first line depending on file version
    if fversion_M == 2
        % Read and parse first line
        line0 = textscan(fgetl(fid),'%2f %*1c %2f %*1c %2f %*1c %2f %*1c %2f %*1c %2f %*1c %5.1f %19.12f %19.12f %19.12f','Whitespace','');

        % Extract SV information
        egnss_i = fgnss_i;
        egnss_f = gnss_i2f(egnss_i);
        PRN = line0{1};
    elseif fversion_M == 3
        % Read and parse first line
        line0 = textscan(fgetl(fid),'%1c %2f %*1c %4f %*1c %2f %*1c %2f %*1c %2f %*1c %2f %*1c %2f %19.12f %19.12f %19.12f','Whitespace','');

        % Extract SV information
        egnss_i = line0{1};
        egnss_f = gnss_i2f(egnss_i);
        PRN = line0{2};
    end
    
    % Select correct array index
    egnss_n = gnss_i2n(egnss_i);

    if any(contains(egnss_i,GNSS_INIT_A))
        % Read and parse broadcast orbit lines
        line1 = textscan(fgetl(fid),'%19.12f',4);
        line2 = textscan(fgetl(fid),'%19.12f',4);
        line3 = textscan(fgetl(fid),'%19.12f',4);
        line4 = textscan(fgetl(fid),'%19.12f',4);
        line5 = textscan(fgetl(fid),'%19.12f',4);
        line6 = textscan(fgetl(fid),'%19.12f',4);
        line7 = textscan(fgetl(fid),'%19.12f',4);

        % Extract remaining SV information
        health = line6{1}(2);
        
        % Extract Keplerian parameters
        M_0 = line1{1}(4); % mean anomaly at reference time [rad]
        e = line2{1}(2); % eccentricity
        sqrt_A = line2{1}(4); % square-root of semi-major axis [m^(1/2)]
        Omega_0 = line3{1}(3); % longitude of ascending node [rad]
        i_0 = line4{1}(1); % inclination angle at reference time [rad]
        omega = line4{1}(3); % argument of perigee [rad]
    
        % Extract Keplerian rate parameters
        Delta_n = line1{1}(3); % mean motion correction [rad/s]
        Omega_dot = line4{1}(4); % rate of right ascension [rad/s]
        i_0_dot = line5{1}(1); % rate of inclination at reference time [rad/s]

        % Extract ephemeris corrections
        C_rc = line4{1}(2); % cosine-harmonic orbit radius correction [m]
        C_rs = line1{1}(2); % sine-harmonic orbit radius correction [m]
        C_ic = line3{1}(2); % cosine-harmonic inclination correction [rad]
        C_is = line3{1}(4); % sine-harmonic inclination correction [rad]
        C_uc = line2{1}(1); % cosine-harmonic argument of latitude correction [rad]
        C_us = line2{1}(3); % sine-harmonic argument of latitude correction [rad]

        % Interpret first line depending on file version
        if fversion_M == 2
            % Extract clock corrections
            a_f0 = line0{8}; % clock bias correction [s]
            a_f1 = line0{9}; % clock drift correction [s/s]
            a_f2 = line0{10}; % clock drift rate correction [s/s^2]

            % Extract reference TOC information (adjustment necessary to account for 2-digit year)
            TOC_Y = line0{2} + floor(TODAY_Y/100)*100 - (line0{2} > mod(TODAY_Y,100))*100; % reference time of clock year [yr]
            TOC = getTimeOfWeek([TOC_Y, line0{3:7}]); % reference time of clock parameters [s]
        elseif fversion_M == 3
            % Extract clock corrections
            a_f0 = line0{9}; % clock bias correction [s]
            a_f1 = line0{10}; % clock drift correction [s/s]
            a_f2 = line0{11}; % clock drift rate correction [s/s^2]

            % Extract reference TOC information
            TOC = getTimeOfWeek([line0{3:8}]); % reference time of clock parameters [s]
        end

        % Extract reference time information
        IODE = line1{1}(1); % issue of data, ephemeris
        TOE = line3{1}(1); % reference time of ephemeris parameters [s]
        TOT = line7{1}(1); % transmission time of message [s]
        WN = line5{1}(3); % week number [wk]

        % Extract other corrections
        URA = line6{1}(1); % user range accuracy (URA) / signal-in-space accuracy (SISA)

        % Append message with generic navigation data
        ephem.(egnss_f)(n_line(egnss_n),1:20) = [PRN, health, ...
            M_0, e, sqrt_A, Omega_0, i_0, omega, ...
            Delta_n, Omega_dot, i_0_dot, ...
            C_rc, C_rs, C_ic, C_is, C_uc, C_us, ...
            a_f0, a_f1, a_f2];
        
        % Extract specific navigation data
        switch egnss_i
            case 'G'
                % Extract reference time information
                IODC = line6{1}(4); % issue of data, clock
                
                % Extract other corrections
                T_GD = line6{1}(3); % L1-L2 group delay correction [s]
                CFI = line7{1}(2); % curve fit interval [h]

                % Extract data/code information
                L2_codes = line5{1}(2); % L2 channel codes
                L2P_flag = line5{1}(4); % L2 P-code data flag

                % Append message with specific navigation data
                ephem.(egnss_f)(n_line(egnss_n),21:N_VARS_PER_GNSS(egnss_n)) = [IODE, IODC, TOE, TOC(1), TOT, WN, ...
                    URA, T_GD, CFI, ...
                    L2_codes, L2P_flag];
            case 'E'
                % Extract reference time information
                IODC = line1{1}(1); % issue of data, clock

                % Extract other corrections
                T_GD_a = line6{1}(3); % E1-E5a group delay correction [s]
                T_GD_b = line6{1}(4); % E1-E5b group delay correction [s]

                % Extract data/code information
                data_sources = line5{1}(2); % data sources

                % Append message with specific navigation data
                ephem.(egnss_f)(n_line(egnss_n),21:N_VARS_PER_GNSS(egnss_n)) = [IODE, IODC, TOE, TOC(1), TOT, WN, ...
                    URA, T_GD_a, T_GD_b, ...
                    data_sources];
            case 'C'
                % Extract reference time information
                IODC = line7{1}(2); % issue of data, clock

                % Extract other corrections
                T_GD_1 = line6{1}(3); % B1I-B3I group delay correction [s]
                T_GD_2 = line6{1}(4); % B2I-B3I group delay correction [s]

                % Append message with specific navigation data
                ephem.(egnss_f)(n_line(egnss_n),21:N_VARS_PER_GNSS(egnss_n)) = [IODE, IODC, TOE, TOC(1), TOT, WN, ...
                    URA, T_GD_1, T_GD_2];
            case 'J'
                % Extract reference time information
                IODC = line6{1}(4); % issue of data, clock

                % Extract other corrections
                T_GD = line6{1}(3); % L1-L2 group delay correction [s]
                if length(line7{:}) == 1
                    CFI = NaN; % curve fit interval (missing)
                else
                    CFI = line7{1}(2); % curve fit interval
                end                

                % Extract data/code information
                L2_codes = line5{1}(2); % L2 channel codes
                L2P_flag = line5{1}(4); % L2 P flag

                % Append message with specific navigation data
                ephem.(egnss_f)(n_line(egnss_n),21:N_VARS_PER_GNSS(egnss_n)) = [IODE, IODC, TOE, TOC(1), TOT, WN, ...
                    URA, T_GD, CFI, ...
                    L2_codes, L2P_flag];
            case 'I'
                % Extract reference time information
                IODC = line1{1}(1); % issue of data, clock

                % Extract other corrections
                T_GD = line6{1}(3); % L1-L2 group delay correction [s]

                % Append message with specific navigation data
                ephem.(egnss_f)(n_line(egnss_n),21:N_VARS_PER_GNSS(egnss_n)) = [IODE, IODC, TOE, TOC(1), TOT, WN, ...
                    URA, T_GD];
        end
    elseif any(contains(egnss_i,GNSS_INIT_B))
        % Read and parse broadcast orbit lines
        line1 = textscan(fgetl(fid),'%19.12f',4);
        line2 = textscan(fgetl(fid),'%19.12f',4);
        line3 = textscan(fgetl(fid),'%19.12f',4);
        if fversion < 3.05
            line4 = {[NaN,NaN,NaN,NaN]};
        else
            line4 = textscan(fgetl(fid),'%19.12f',4);
        end

        % Extract remaining SV information
        health = line1{1}(4);

        % Extract satellite position information (PZ-90 coordinate system)
        sat_x = line1{1}(1)*1e3; % satellite x-axis position [m]
        sat_y = line2{1}(1)*1e3; % satellite y-axis position [m]
        sat_z = line3{1}(1)*1e3; % satellite z-axis position [m]

        % Extract satellite velocity information
        sat_u = line1{1}(2)*1e3; % satellite x-axis velocity [m/s]
        sat_v = line2{1}(2)*1e3; % satellite y-axis velocity [m/s]
        sat_w = line3{1}(2)*1e3; % satellite z-axis velocity [m/s]

        % Extract satellite acceleration information
        sat_u_dot = line1{1}(3)*1e3; % satellite x-axis acceleration [m/s^2]
        sat_v_dot = line2{1}(3)*1e3; % satellite y-axis acceleration [m/s^2]
        sat_w_dot = line3{1}(3)*1e3; % satellite z-axis acceleration [m/s^2]

        % Append message with generic navigation data
        ephem.(egnss_f)(n_line(egnss_n),1:11) = [PRN, health, ...
            sat_x, sat_y, sat_z, sat_u, sat_v, sat_w, sat_u_dot, sat_v_dot, sat_w_dot];

        % Extract specific navigation data
        switch egnss_i
            case 'R'
                % Extract remaining SV information
                status_flags = line4{1}(1); % status flags (9-bit register)
                health_flags = line4{1}(4); % health flags (9-bit register)

                % Extract clock corrections
                tau_n = line0{9}; % clock bias correction [s]
                gamma_n = line0{10}; % relative frequency bias correction [s/s]

                % Extract reference time information
                t_k = line0{11}; % message frame time [s]
                AOII = line3{1}(4); % age of immediate information [day]

                % Extract frequency information
                freq_num = line2{1}(4); % frequency (channel) number of transmitted signals

                % Extract other corrections
                T_GD = line4{1}(2); % L1-L2 group delay correction [s]
                URA = line4{1}(3); % user range accuracy (URA)

                % Append message with specific navigation data
                ephem.(egnss_f)(n_line(egnss_n),12:N_VARS_PER_GNSS(egnss_n)) = [tau_n, gamma_n, ...
                    t_k, AOII, freq_num, ...
                    T_GD, URA, ...
                    status_flags, health_flags];
            case 'S'
                % Extract clock corrections
                a_f0 = line0{9}; % clock bias correction [s]
                a_f1 = line0{10}; % relative frequency bias correction [s/s]

                % Extract reference time information
                IODE = line3{1}(4); % issue of data, ephemeris
                IODC = line3{1}(4); % issue of data, clock
                TOT = line0{11}; % transmission time of message [s]
                WN = NaN; % week number [wk]

                % Extract other corrections
                URA = line2{1}(4); % user range accuracy (URA)

                % Append message with specific navigation data
                ephem.(egnss_f)(n_line(egnss_n),12:N_VARS_PER_GNSS(egnss_n)) = [a_f0, a_f1, ...
                    IODE, IODC, TOT, WN, ...
                    URA];
        end
    end

    n_line(egnss_n) = n_line(egnss_n) + 1;
end

% Close file
fclose(fid);

% Remove outdated entries following off-hour cut-over if requested
if opts.CleanOutput
    for s = 1:N_GNSS
        if N_entries_per_gnss(s) ~= 0
            PRN_unique = unique(ephem.(GNSS_FULL{s})(:,1))';
            for egnss_n = PRN_unique
                PRN_indices = find(ephem.(GNSS_FULL{s})(:,1) == egnss_n);
                if ~isempty(PRN_indices)
                    new_uploads = find(mod(ephem.(GNSS_FULL{s})(PRN_indices,23),3600) ~= 0);
                    new_uploads(new_uploads == PRN_indices(end)) = [];
                    if ~isempty(new_uploads)
                        ephem.(GNSS_FULL{s})((ephem.(GNSS_FULL{s})(new_uploads + 1,23) - ephem.(GNSS_FULL{s})(new_uploads,23)) < 600,:) = [];
                    end
                end
            end
        end
    end
end

% Convert array to struct/table if requested
switch opts.OutputType
    case 'struct'
        for s = 1:N_GNSS
            if N_entries_per_gnss(s) ~= 0
                ephem_temp = ephem.(GNSS_FULL{s});
                ephem.(GNSS_FULL{s}) = struct;
                for v = 1:N_VARS_PER_GNSS(s)
                    ephem.(GNSS_FULL{s}).(VARIABLE_NAMES.(GNSS_FULL{s})(v)) = ephem_temp(:,v);
                end
            end
        end
    case 'table'
        for s = 1:N_GNSS
            if N_entries_per_gnss(s) ~= 0
                ephem_temp = ephem.(GNSS_FULL{s});
                ephem.(GNSS_FULL{s}) = table;
                for v = 1:N_VARS_PER_GNSS(s)
                    ephem.(GNSS_FULL{s}).(VARIABLE_NAMES.(GNSS_FULL{s})(v)) = ephem_temp(:,v);
                end
                ephem.(GNSS_FULL{s}).Properties.VariableUnits = VARIABLE_UNITS.(GNSS_FULL{s});
            end
        end
end

% Pull data to top-level (replacing top-level struct) if single-GNSS file
if fgnss_i ~= 'M'
    ephem = ephem.(fgnss_f);
    if contains(fgnss_f,fieldnames(iono))
        iono = iono.(fgnss_f);
    end
    if contains(fgnss_f,fieldnames(time))
        time = time.(fgnss_f);
    end
end

end
