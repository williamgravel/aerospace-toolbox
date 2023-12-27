function [params] = TRF(gnss)
arguments
    gnss (1,1) string {mustBeValidSatelliteSystem} = 'gps'
end

switch gnss
    case 'gps'
        model = 'WGS 84';
        a = 6378137.0; % equatorial radius [m]
        f = 1/298.257223563; % flattening
        Omega_dot = 7.2921151467e-5; % rotation rate [rad/s]
        mu = 3.986005e14; % gravitational constant [m^3/s^2]
        c = 2.99792458e8; % speed of light [m/s]
        J_2 = 0.0010826262; % second zonal harmonic coefficient
    case 'glonass'
        model = 'PZ-90.02';
        a = 6378136; % equatorial radius [m]
        f = 1/298.25784; % flattening
        Omega_dot = 7.292115e-5; % rotation rate [rad/s]
        mu = 3.986004418e14; % gravitational constant [m^3/s^2]
        c = 2.99792458e8; % speed of light [m/s]
        J_2 = 0.00108262575; % second zonal harmonic coefficient
    case 'galileo'
        model = 'GTRF19v01';
        a = 6378137.0; % equatorial radius [m]
        f = 1/298.257222101; % flattening
        Omega_dot = 7.2921151467e-5; % rotation rate [rad/s]
        mu = 3.986004418e14; % gravitational constant [m^3/s^2]
        c = 2.99792458e8; % speed of light [m/s]
        J_2 = NaN; % second zonal harmonic coefficient
    case 'beidou'
        model = 'CTRF2000';
        a = 6378137.0; % equatorial radius [m]
        f = 1/298.257222101; % flattening
        Omega_dot = 7.292115e-5; % rotation rate [rad/s]
        mu = 3.986004418e14; % gravitational constant [m^3/s^2]
        c = 2.99792458e8; % speed of light [m/s]
        J_2 = 0.00108263; % second zonal harmonic coefficient
end

e = sqrt(2*f - f^2); % eccentricity
b = a*(1 - f); % polar radius [m]

params = struct('model',model,'a',a,'f',f,'e',e,'b',b,'mu',mu,'Omega_dot',Omega_dot,'c',c,'J_2',J_2);

end
