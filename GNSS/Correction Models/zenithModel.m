function [ZTD,ZHD,ZWD,a_h,a_w] = zenithModel(lat,lon,h,doy,opts)
arguments
    lat (1,1) double
    lon (1,1) double
    h (1,1) double
    doy (1,1) double
    opts.MeteoModel (1,1) string
    opts.DelayModel (1,1) string
    opts.RefinedDecimals (1,1) logical
end

% Define day of year angle
doy_p = (doy - 28)/365.25*2*pi; % day of year angle [rad]

% Define molar properties of dry air and water
M_d = 28.9644e-3; % molar mass of dry air [kg/mol]
M_w = 18.0152e-3; % molar mass of water [kg/mol]
R = 8.314; % molar gas constant [J/(mol*K)]
R_d = R/M_d; % gas constant of dry air [J/(kg*K)]

% Define refractivity constants
k_1 = 77.604; % refractivity constant 1 [K/mbar]
k_2 = 64.79; % refractivity constant 2 [K/mbar]
k_3 = 3.776e5; % refractivity constant 3 [K^2/mbar]
k_2_p = k_2 - k_1*M_w/M_d; % modified refractivity constant 2 [K/mbar]

% Define gravitational potential function
f_g = @(lat,h) 1 - 0.00266*cosd(2*lat) - 0.00000028*h; % where lat = geodetic latitude [deg], h = ellipsoidal height [m]

% Calculate gravitational acceleration at altitude
g_m_0 = 9.80665; % mean gravitational acceleration at MSL [m/s^2]
g_m = g_m_0/f_g(lat,h); % mean gravitational acceleration [m/s^2]

switch opts.MeteoModel
    case "SZTD"
        SZTD = load("Lookup Tables\SZTD.mat");
        
        P_0 = interp1(SZTD.phi,SZTD.a_0.P_0,lat) + interp1(SZTD.phi,SZTD.A.P_0,lat)*cos(doy_p); % air pressure [mbar]
        T_0 = interp1(SZTD.phi,SZTD.a_0.T_0,lat) + interp1(SZTD.phi,SZTD.A.T_0,lat)*cos(doy_p); % temperature [K]
        e_0 = interp1(SZTD.phi,SZTD.a_0.e_0,lat) + interp1(SZTD.phi,SZTD.A.e_0,lat)*cos(doy_p); % partial pressure of water vapor [mbar]
        beta = interp1(SZTD.phi,SZTD.a_0.beta,lat) + interp1(SZTD.phi,SZTD.A.beta,lat)*cos(doy_p); % temperature lapse rate [K/m]
        lambda = interp1(SZTD.phi,SZTD.a_0.lambda,lat) + interp1(SZTD.phi,SZTD.A.lambda,lat)*cos(doy_p); % water vapor decrease factor
        T_m = interp1(SZTD.phi,SZTD.a_0.T_m,lat) + interp1(SZTD.phi,SZTD.A.T_m,lat)*cos(doy_p); % mean temperature of water vapor [K]

        a_h = NaN;
        a_w = NaN;
    case "UNB3m" % only valid for latitudes between 15 and 75 degrees
        UNB3m = load("Lookup Tables\UNB3m.mat");
        
        P_0 = interp1(UNB3m.phi,UNB3m.a_0.P_0,lat) + interp1(UNB3m.phi,UNB3m.A.P_0,lat)*cos(doy_p); % air pressure [mbar]
        T_0 = interp1(UNB3m.phi,UNB3m.a_0.T_0,lat) + interp1(UNB3m.phi,UNB3m.A.T_0,lat)*cos(doy_p); % temperature [K]
        rh = interp1(UNB3m.phi,UNB3m.a_0.rh,lat) + interp1(UNB3m.phi,UNB3m.A.rh,lat)*cos(doy_p); % relative humidity [%]
        beta = interp1(UNB3m.phi,UNB3m.a_0.beta,lat) + interp1(UNB3m.phi,UNB3m.A.beta,lat)*cos(doy_p); % temperature lapse rate [K/m]
        lambda = interp1(UNB3m.phi,UNB3m.a_0.lambda,lat) + interp1(UNB3m.phi,UNB3m.A.lambda,lat)*cos(doy_p); % water vapor decrease factor

        e_s = 0.01*exp(1.2378847e-5*T_0^2 - 1.9121316e-2*T_0 + 33.93711047 - 6.3431645e3*T_0^(-1)); % saturation pressure of water vapor [mbar]
        f_w = 1.00062 + 3.14e-6*P_0 + 5.6e-7*(T_0 - 273.15)^2; % enhancement factor
        e_0 = rh/100*e_s*f_w; % partial pressure of water vapor [mbar]

        T_m = (T_0 - beta*h)*(1 - beta*R_d/(g_m*(lambda + 1)));

        a_h = NaN;
        a_w = NaN;
    case "GPT3"
        GPT3 = gpt3_5_fast_readGrid;

        [P_0,T_0,beta,T_m,e_0,a_h,a_w,lambda] = gpt3_5_fast(deg2rad(lat),deg2rad(lon),h,doy,0,GPT3);
        T_0 = T_0 + 273.15;
end

% Define pressure ratio to temperature ratio conversion exponent
mu = g_m/(R_d*beta) - 1; % thermodynamic ratio exponent
mu_p = mu + 1; % modified thermodynamic ratio exponent

% Define modified water vapor decrease factor
lambda_p = lambda + 1; % modified water vapor decrease factor

% Define temperature at altitude (only valid in first layer of troposphere, under 11 km)
T = T_0 - beta*h; % temperature [K]

% Compute Zenith delays
ZHD = 0.0022768/f_g(lat,h)*P_0*(T/T_0)^mu_p; % Zenith hydrostatic delay [m]
ZWD = 1e-6*(k_2_p + k_3/T_m)*R_d/(lambda_p*g_m)*e_0*(T/T_0)^(mu_p*lambda_p); % Zenith wet delay [m]
ZTD = ZHD + ZWD; % Zenith total delay [m]

end
