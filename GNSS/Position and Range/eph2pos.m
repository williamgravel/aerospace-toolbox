function [pos,health,b,epsilon_r] = eph2pos(ephem,t,PRN,opts)
arguments
    ephem (1,1) struct
    t (:,2) double
    PRN (1,1) double
    opts.ReferenceFrame (1,1) string = 'gps'
end

% Fetch defining contants for ellipsoid model
model = TRF(opts.ReferenceFrame);

% Initialize empty arrays
sz = size(t,1);
E_k = zeros(sz,1);
pos = NaN(sz,3);
health = NaN(sz,1);
b = NaN(sz,1);
epsilon_r = NaN(sz,1);

% Determine target PRN indices
ii = find(ephem.PRN == PRN);

if isempty(ii)
    warning('Unable to find ephemeris entries with requested PRN. Terminating program.')
    return
end

% Define week to seconds conversion factor
W2S = 7*24*60*60; % week to seconds conversion factor [s/wk]

% Compute elapsed time of each TOE/TOC w.r.t. first ephemeris entry
TOE_norm = (ephem.WN(ii) - ephem.WN(ii(1)))*W2S + (ephem.TOE(ii) - ephem.TOE(ii(1)));
TOC_norm = (ephem.WN(ii) - ephem.WN(ii(1)))*W2S + (ephem.TOC(ii) - ephem.TOC(ii(1)));

% Compute elapsed times of each input time w.r.t. first ephemeris entry
dTOE_norm = (t(:,2) - ephem.WN(ii(1)))*W2S + (t(:,1) - ephem.TOE(ii(1)));
dTOC_norm = (t(:,2) - ephem.WN(ii(1)))*W2S + (t(:,1) - ephem.TOC(ii(1)));

% Find indices of nearest ephemeris entries in time to input time
[~,jj] = min(abs(repmat(dTOE_norm,1,size(TOE_norm,1)) - TOE_norm'),[],2);

% Combine PRN- and time-matching indices
kk = ii(jj);

% Determine time from ephemeris and clock reference epochs [s]
t_ke = dTOE_norm - TOE_norm(jj);
t_kc = dTOC_norm - TOC_norm(jj);

% Compute corrected mean motion [rad/s]
n_0 = sqrt(model.mu./ephem.sqrt_A(kk).^6);
n_A = n_0 + ephem.Delta_n(kk);

% Compute mean anomaly [rad]
Delta_M = n_A.*t_ke;
M_k = ephem.M_0(kk) + Delta_M;

% Compute eccentric anomaly [rad]
for p = 1:sz
    E_k(p) = mean2eccentric(M_k(p),ephem.e(kk(p)),'AngleUnits','rad');
end

% Compute true anomaly [rad]
nu_k = 2*atan2(sqrt(1 + ephem.e(kk)).*tan(E_k/2),sqrt(1 - ephem.e(kk)));

% Compute argument of latitude [rad]
Phi_k = nu_k + ephem.omega(kk);

% Determine second harmonic perturbations for argument of latitude [rad], radius [m], and inclination [rad]
delta_u_k = ephem.Cus(kk).*sin(2*Phi_k) + ephem.Cuc(kk).*cos(2*Phi_k);
delta_r_k = ephem.Crs(kk).*sin(2*Phi_k) + ephem.Crc(kk).*cos(2*Phi_k);
delta_i_k = ephem.Cis(kk).*sin(2*Phi_k) + ephem.Cic(kk).*cos(2*Phi_k);

% Compute corrected argument of latitude [rad], radius [m], and inclination [rad]
u_k = Phi_k + delta_u_k;
r_k = ephem.sqrt_A(kk).^2.*(1 - ephem.e(kk).*cos(E_k)) + delta_r_k;
i_k = ephem.i_0(kk) + ephem.i_0_dot(kk).*t_ke + delta_i_k;

% Compute positions in orbital plane [m]
x_k_p = r_k.*cos(u_k);
y_k_p = r_k.*sin(u_k);

% Compute corrected longitude of ascending node [rad]
Omega_k = ephem.Omega_0(kk) + (ephem.Omega_dot(kk) - model.Omega_dot).*t_ke - model.Omega_dot*ephem.TOE(kk);

% Compute positions in Earth-fixed coordinates [m]
x_k = x_k_p.*cos(Omega_k) - y_k_p.*cos(i_k).*sin(Omega_k);
y_k = x_k_p.*sin(Omega_k) + y_k_p.*cos(i_k).*cos(Omega_k);
z_k = y_k_p.*sin(i_k);

% Assemble SV positions and health
pos = [x_k, y_k, z_k];
health = ephem.health(kk);

% Compute SV clock bias and relativity corrections [m]
b = model.c*(ephem.a_f0(kk) + ephem.a_f1(kk).*t_kc + ephem.a_f2(kk).*t_kc.^2);
F = -2*sqrt(model.mu)/model.c^2;
epsilon_r = model.c*(F*ephem.e(kk).*ephem.sqrt_A(kk).*sin(E_k));

end
