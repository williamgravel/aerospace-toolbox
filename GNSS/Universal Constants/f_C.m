function [f] = f_C(gnss,signal,channel)
arguments
    gnss (1,1) string {mustBeValidSatelliteSystem}
    signal (1,1) string
    channel (1,1) double = 0
end

switch gnss
    case 'gps'
        switch signal
            case 'L1'
                f = 1575.42;
            case 'L2'
                f = 1227.60;
            case 'L5'
                f = 1176.45;
            otherwise
                f = NaN;
        end
    case 'glonass'
        switch signal
            case 'L1'
                f = 1602 + channel*0.5625;
            case 'L2'
                f = 1246 + channel*0.4375;
            case 'L3'
                f = 1202.025;
            otherwise
                f = NaN;
        end
    case 'galileo'
        switch signal
            case 'E1'
                f = 1575.420;
            case 'E5'
                f = 1191.795;
            case 'E5a'
                f = 1176.450;
            case 'E5b'
                f = 1207.140;
            case 'E6'
                f = 1278.750;
            otherwise
                f = NaN;
        end
    case 'beidou'
        switch signal
            case 'B1I'
                f = 1561.098;
            case 'B1Q'
                f = 1561.098;
            case 'B1C'
                f = 1575.420;
            case 'B1A'
                f = 1575.420;
            case 'B2I'
                f = 1207.140;
            case 'B2Q'
                f = 1207.140;
            case 'B2a'
                f = 1176.450;
            case 'B2b'
                f = 1207.140;
            case 'B3I'
                f = 1268.520;
            case 'B3Q'
                f = 1268.520;
            case 'B3A'
                f = 1268.520;
        end
    case 'qzss'
        switch signal
            case 'L1'
                f = 1575.42;
            case 'L2'
                f = 1227.60;
            case 'L5'
                f = 1176.45;
            case 'L6'
                f = 1278.75;
            otherwise
                f = NaN;
        end
    case 'irnss'
        switch signal
            case 'L5'
                f = 1176.450;
            case 'S'
                f = 2492.028;
            otherwise
                f = NaN;
        end
end

f = f*1e6;

end
