function [ephem,iono] = readAlmanac(filename,opts)
arguments
    filename (1,1) string
    opts.OutputType (1,1) string {mustBeMember(opts.OutputType,{'struct','array'})} = 'struct'
    opts.WeekRollovers (1,1) double = 0
end

fid = fopen(filename,'r');

alpha_0 = 0; % ionospheric vertical delay amplitude coefficient [s]
alpha_1 = 0; % ionospheric vertical delay amplitude coefficient [s/semi-circle]
alpha_2 = 0; % ionospheric vertical delay amplitude coefficient [s/(semi-circle)^2]
alpha_3 = 0; % ionospheric vertical delay amplitude coefficient [s/(semi-circle)^3]
beta_0 = 0; % ionospheric period coefficient [s]
beta_1 = 0; % ionospheric period coefficient [s/semi-circle]
beta_2 = 0; % ionospheric period coefficient [s/(semi-circle)^2]
beta_3 = 0; % ionospheric period coefficient [s/(semi-circle)^3]

fid_count = fopen(filename,'r');
N_entries = length(strfind(fread(fid_count)',[80,82,78]));
fclose(fid_count);

ephem_array = zeros(N_entries,29);
i = 1;

% Data record
while true
    if feof(fid)
        break
    end

    % Read and parse lines
    fgetl(fid);
    line1 = textscan(fgetl(fid),'%*25s %f','Delimiter','');
    line2 = textscan(fgetl(fid),'%*25s %f','Delimiter','');
    line3 = textscan(fgetl(fid),'%*25s %f','Delimiter','');
    line4 = textscan(fgetl(fid),'%*25s %f','Delimiter','');
    line5 = textscan(fgetl(fid),'%*25s %f','Delimiter','');
    line6 = textscan(fgetl(fid),'%*25s %f','Delimiter','');
    line7 = textscan(fgetl(fid),'%*25s %f','Delimiter','');
    line8 = textscan(fgetl(fid),'%*25s %f','Delimiter','');
    line9 = textscan(fgetl(fid),'%*25s %f','Delimiter','');
    line10 = textscan(fgetl(fid),'%*25s %f','Delimiter','');
    line11 = textscan(fgetl(fid),'%*25s %f','Delimiter','');
    line12 = textscan(fgetl(fid),'%*25s %f','Delimiter','');
    line13 = textscan(fgetl(fid),'%*25s %f','Delimiter','');
    fgetl(fid);

    % Extract SV information
    PRN = line1{1};
    health = line2{1};

    % Extract time information
    TOC = 0; % reference time of clock parameters [s]
    TOE = line4{1}; % reference time of ephemeris parameters [s]
    IODE = 0; % issue of data, ephemeris
    IODC = 0; % issue of data, clock
    TOT = 0; % transmission time of message [s]
    WN = line13{1} + opts.WeekRollovers*1024; % GPS week number

    % Extract Keplerian parameters
    M_0 = line10{1}; % mean anomaly at reference time [rad]
    e = line3{1}; % eccentricity
    sqrt_A = line7{1}; % square-root of semi-major axis [m^(1/2)]
    Omega_0 = line8{1}; % longitude of ascending node at reference time [rad]
    i_0 = line5{1}; % inclination angle at reference time [rad]
    omega = line9{1}; % argument of perigee [rad]

    % Extract Keplerian rate parameters
    Delta_n = 0; % mean motion correction [rad/s]
    Omega_dot = line6{1}; % rate of right ascension [rad/s]
    i_0_dot = 0; % rate of inclination at reference time [rad/s]

    % Extract ephemeris corrections
    Crc = 0; % cosine-harmonic orbit radius correction [m]
    Crs = 0; % sine-harmonic orbit radius correction [m]
    Cic = 0; % cosine-harmonic inclination correction [rad]
    Cis = 0; % sine-harmonic inclination correction [rad]
    Cuc = 0; % cosine-harmonic argument of latitude correction [rad]
    Cus = 0; % sine-harmonic argument of latitude correction [rad]

    % Extract clock corrections
    a_f0 = line11{1}; % clock bias correction [s]
    a_f1 = line12{1}; % clock drift correction [s/s]
    a_f2 = 0; % clock drift rate correction [s/s^2]

    % Extract other corrections
    TGD = 0; % group delay correction [s]
    CFI = 0; % curve fit interval [h]
    URA_N = 0; % user range accuracy index
    if URA_N < 6
        URA = round(2^(1 + URA_N/2),1); % user range accuracy [m]
    elseif URA_N < 15
        URA = 2^(URA_N - 2); % user range accuracy [m]
    else
        URA = 0; % user range accuracy [m]
    end

    % Append message
    ephem_array(i,:) = [PRN, M_0, Delta_n, e, sqrt_A, Omega_0, i_0, omega, Omega_dot, i_0_dot, ...
        Cuc, Cus, Crc, Crs, Cic, Cis, TOE, IODE, WN, TOC, a_f0, a_f1, a_f2, TGD, health, ...
        IODC, TOT, CFI, URA];
    i = i + 1;
end

fclose(fid);

if strcmpi(opts.OutputType,'array')
    ephem = ephem_array(:,1:25);
    iono = [alpha_0, alpha_1, alpha_2, alpha_3, beta_0, beta_1, beta_2, beta_3];
elseif strcmpi(opts.OutputType,'struct')
    ephem = struct('PRN',ephem_array(:,1),'health',ephem_array(:,25),'TOC',ephem_array(:,20),'TOE',ephem_array(:,17),...
        'IODE',ephem_array(:,18),'IODC',ephem_array(:,26),'t_0t',ephem_array(:,27),'WN',ephem_array(:,19),...
        'M_0',ephem_array(:,2),'e',ephem_array(:,4),'sqrt_A',ephem_array(:,5),'Omega_0',ephem_array(:,6),...
        'i_0',ephem_array(:,7),'omega',ephem_array(:,8),'Delta_n',ephem_array(:,3),'Omega_dot',ephem_array(:,9),...
        'i_0_dot',ephem_array(:,10),'Crc',ephem_array(:,13),'Crs',ephem_array(:,14),'Cic',ephem_array(:,15),...
        'Cis',ephem_array(:,16),'Cuc',ephem_array(:,11),'Cus',ephem_array(:,12),'a_f0',ephem_array(:,21),...
        'a_f1',ephem_array(:,22),'a_f2',ephem_array(:,23),'T_GD',ephem_array(:,24),'CFI',ephem_array(:,28),...
        'URA',ephem_array(:,29));
    iono = struct('alpha_0',alpha_0,'alpha_1',alpha_1,'alpha_2',alpha_2,'alpha_3',alpha_3,...
        'beta_0',beta_0,'beta_1',beta_1,'beta_2',beta_2,'beta_3',beta_3);
end

end
