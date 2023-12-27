function [obs,lli,rss,header] = readObservation(filename,opts)
arguments
    filename (1,1) string
    opts.OutputType (1,1) string {mustBeMember(opts.OutputType,{'struct','array'})} = 'struct'
    opts.SplitConstellations (1,1) logical = false
end

fid = fopen(filename,'r');

cline = '';

header = struct;

% Header section
while ~contains(cline,'END OF HEADER')
    cline = fgetl(fid);
    clabel = cline(61:end);

    if contains(cline(61:end),'COMMENT')
        continue
    elseif contains(clabel,'RINEX VERSION / TYPE')
        hline1 = textscan(cline,'%9.2f %*11c %1c %*19c %1c %*19c','Delimiter','\n','Whitespace','');
    elseif contains(clabel,'PGM / RUN BY / DATE')

    elseif contains(clabel,'MARKER NAME')

    elseif contains(clabel,'OBSERVER / AGENCY')

    elseif contains(clabel,'REC # / TYPE / VERS')

    elseif contains(clabel,'ANT # / TYPE')

    elseif contains(clabel,'APPROX POSITION XYZ')

    elseif contains(clabel,'ANTENNA: DELTA H/E/N')

    elseif contains(clabel,'# / TYPES OF OBSERV')
        N_obs = str2double(cline(1:6));
        T_obs = [];
        for i = 0:floor(N_obs/9)
            if ~isempty(T_obs)
                cline = fgetl(fid);
            end
            T_obs_cell = textscan(cline(7:60),'%*4s %2s',9,'Delimiter','\n','Whitespace','','TextType','string');
            T_obs = [T_obs; T_obs_cell{:}];
        end
        T_obs = T_obs(1:N_obs)';
    elseif contains(clabel,'# OF SATELLITES')

    elseif contains(clabel,'PRN / # OF OBS')

    elseif contains(clabel,'INTERVAL')

    elseif contains(clabel,'TIME OF FIRST OBS')

    elseif contains(clabel,'TIME OF LAST OBS')

    end
end

% O_codes = ["C1","C2","C5","C6","C7","C8","L1","L2","L5","L6","L7","L8","D1","D2","D5","D6","D7","D8","S1","S2","S5","S6","S7","S8","P1","P2"];
% O_lines = floor((N_obs - 1)/5);

obs_array = zeros(0,N_obs);
lli_array = zeros(0,N_obs,'int8');
rss_array = zeros(0,N_obs,'int8');
prn_array = zeros(0,1);
epoch_array = zeros(0,2);

L_per_SV = ceil(N_obs/5);

i = 1;

% Data record
while true
    if feof(fid)
        break
    end

    cline = fgetl(fid);
    lineA = textscan(cline(1:32),'%3f %*1s %2f %*1s %2f %*1s %2f %*1s %2f %11.7f %*2s %1f %3f','Delimiter','\n','Whitespace','');
    lineB = textscan(cline(33:end),'%3s',12,'TextType','string');
    sv = lineB{:};

    Y = 2000 + lineA{1};
    M = lineA{2};
    D = lineA{3};
    H = lineA{4};
    MI = lineA{5};
    S = floor(lineA{6});
    [tow,week] = formatGNSSTime('GPS',Y,M,D,H,MI,S);
    epoch_flag = lineA{7};
    SV_per_E = lineA{8};

    for p = 1:(ceil(SV_per_E/12) - 1)
        cline = fgetl(fid);
        lineB = textscan(cline(33:end),'%3s',12,'TextType','char');
        sv = [sv; lineB{:}];
    end

    obs_array = [obs_array; zeros(SV_per_E,N_obs)];
    lli_array = [lli_array; zeros(SV_per_E,N_obs,'int8')];
    rss_array = [rss_array; zeros(SV_per_E,N_obs,'int8')];
    prn_array = [prn_array; zeros(SV_per_E,1)];
    epoch_array = [epoch_array; ones(SV_per_E,2).*[week,tow]];

    % Read and parse lines
    for s = 1:SV_per_E
        switch sv{s}(1)
            case 'G'
                prn_array(i) = str2double(sv{s}(2:3));
            case 'R'
                prn_array(i) = str2double(sv{s}(2:3)) + 100;
            case 'S'
                prn_array(i) = str2double(sv{s}(2:3)) + 200;
            case 'E'
                prn_array(i) = str2double(sv{s}(2:3)) + 300;
            case 'C'
                prn_array(i) = str2double(sv{s}(2:3)) + 400;
            case 'J'
                prn_array(i) = str2double(sv{s}(2:3)) + 500;
            case 'I'
                prn_array(i) = str2double(sv{s}(2:3)) + 600;
            otherwise
                prn_array(i) = str2double(sv{s}(2:3)) + 700;
        end

        for c = 1:L_per_SV
            cline = fgetl(fid);
            cline(cline == ' ') = '0';
            lineC = textscan(cline,'%14.3f %1u8 %1u8',5);
            O_num = (1:5) + 5*(c - 1);
            obs_array(i,O_num(O_num <= N_obs)) = lineC{1}';
            lli_array(i,O_num(O_num <= N_obs)) = lineC{2}';
            rss_array(i,O_num(O_num <= N_obs)) = lineC{3}';
        end

        i = i + 1;
    end
end

fclose(fid);

if strcmpi(opts.OutputType,'array')
    obs = [prn_array, epoch_array, obs_array];
    lli = [prn_array, epoch_array, lli_array];
    rss = [prn_array, epoch_array, rss_array];
elseif strcmpi(opts.OutputType,'struct')
    if opts.SplitConstellations

    else

    end
    
    obs = struct('PRN',prn_array,'week',epoch_array(:,1),'tow',epoch_array(:,2));
    lli = struct('PRN',prn_array,'week',epoch_array(:,1),'tow',epoch_array(:,2));
    rss = struct('PRN',prn_array,'week',epoch_array(:,1),'tow',epoch_array(:,2));

    for i = 1:length(T_obs)
        obs.(T_obs(i)) = obs_array(:,i);
        lli.(T_obs(i)) = lli_array(:,i);
        rss.(T_obs(i)) = rss_array(:,i);
    end
end

end
