function [n_N2,n_O,n_O2,n_Ar,n_He,n_H] = CalcNumberDensity(Z,T,u,t)
%CALCNUMBERDENSITY Calculates gas-specific number densities for altitudes above boundary layer.
%
%   Inputs:     Z                       -   geometric altitude
%               T                       -   kinetic temperature
%               u                       -   universal parameters structure
%               t                       -   tables structure
%
%   Outputs:    n_N2                    -   molecular nitrogen number density
%               n_O                     -   atomic oxygen number density
%               n_O2                    -   molecular oxygen number density
%               n_Ar                    -   atomic argon number density
%               n_He                    -   atomic helium number density
%               n_H                     -   atomic hydrogen number density
%
%   Author:     William Gravel
%   Created:    05/06/2021
%   Edited:     05/19/2021
%   Purpose:    COESA's U.S. Standard Atmosphere 1976 Model

%% Setup
% Extract universal constants, equations, and tables
c = u.c;
e = u.e;
Atmos = t.Atmos;
AirComp = t.AirComp;
DiffCoeffs = t.DiffCoeffs;
FluxCoeffs = t.FluxCoeffs;

%% Reference Values
% Define reference altitude vector
Z_ref_1 = 86:0.005:90;
Z_ref_2 = 90:0.01:300;
Z_ref_3 = 300:0.02:500;
Z_ref_4 = 500:0.05:1000;
if max(Z) < 500
    Z_ref = [Z_ref_1'; Z_ref_2(2:end)'; Z_ref_3(2:end)'];
else
    Z_ref = [Z_ref_1'; Z_ref_2(2:end)'; Z_ref_3(2:end)'; Z_ref_4(2:end)'];
end

% Compute height-dependent quantities for reference altitudes
g_ref = e.g(Z_ref)*1000;
lambda = Atmos.Z{'9','L_Kb'}/(c.T_inf - Atmos.Z{'10','T_b'});
zeta_ref = (Z_ref - Atmos.Z{'10','Z_b'})*(c.r_0 + Atmos.Z{'10','Z_b'})./(c.r_0 + Z_ref);
T_ref = T_Z(Z_ref);
dT_dZ_ref = dT_dZ_Z(Z_ref);
K_ref = K_Z(Z_ref);

% Define first mean molecular weight approximation
M_1 = zeros(length(Z_ref),1);
M_1(Z_ref <= 100) = c.M_0;
M_1(Z_ref > 100) = AirComp.BL{'N2','M_i'};

% Define species-dependent number density prefix function
n = @(G) AirComp.BL{G,'n_i'}*Atmos.Z{'7','T_b'}./T_ref;

%% Molecular Nitrogen
% Calculate molecular nitrogen number density for reference altitudes
tau_N2_ref = cumtrapz(Z_ref,M_1.*g_ref./(c.R_star*T_ref));
n_N2_ref = n('N2').*exp(-tau_N2_ref);

%% First Background Gas Approximation
% Determine first total number density approximation
sum_n_1 = n_N2_ref;

%% Atomic Oxygen
% Calculate integral term for atomic oxygen number density
D_O = D_i('O',sum_n_1);
f_O = f_i('O',D_O,M_1);
v_O = v_i('O',Z_ref);
tau_O_ref = cumtrapz(Z_ref,f_O + v_O);
n_O_ref = n('O').*exp(-tau_O_ref);

%% Molecular Oxygen
% Calculate integral term for molecular oxygen number density
D_O2 = D_i('O2',sum_n_1);
f_O2 = f_i('O2',D_O2,M_1);
v_O2 = v_i('O2',Z_ref);
tau_O2_ref = cumtrapz(Z_ref,f_O2 + v_O2);
n_O2_ref = n('O2').*exp(-tau_O2_ref);

%% Second Background Gas Approximation
% Define second background gas approximation
sum_n_2 = n_N2_ref + n_O_ref + n_O2_ref;
M_2 = zeros(length(Z_ref),1);
M_2(Z_ref <= 100) = c.M_0;
if any(Z_ref > 100)
    M_2(Z_ref > 100) = sum([n_N2_ref(Z_ref > 100),n_O_ref(Z_ref > 100),n_O2_ref(Z_ref > 100)].*AirComp.BL{{'N2','O','O2'},'M_i'}',2)./sum_n_2(Z_ref > 100);
end

%% Atomic Argon
% Calculate integral term for atomic argon number density
D_Ar = D_i('Ar',sum_n_2);
f_Ar = f_i('Ar',D_Ar,M_2);
v_Ar = v_i('Ar',Z_ref);
tau_Ar_ref = cumtrapz(Z_ref,f_Ar + v_Ar);
n_Ar_ref = n('Ar').*exp(-tau_Ar_ref);

%% Atomic Helium
% Calculate integral term for atomic helium number density
D_He = D_i('He',sum_n_2);
f_He = f_i('He',D_He,M_2);
v_He = v_i('He',Z_ref);
tau_He_ref = cumtrapz(Z_ref,f_He + v_He);
n_He_ref = n('He').*exp(-tau_He_ref);

%% Third Background Gas Approximation
% Define third total number density approximation
sum_n_3 = n_N2_ref + n_O_ref + n_O2_ref + n_Ar_ref + n_He_ref;

%% Atomic Hydrogen
% Calculate atomic hydrogen number density
i_Z_11 = find(Z_ref == Atmos.Z{'11','Z_b'},1);
i_Z_150 = find(Z_ref == 150,1);

D_H = DiffCoeffs{'H','a_i'}./sum_n_3(i_Z_150:end).*(T_ref(i_Z_150:end)/273.15).^DiffCoeffs{'H','b_i'};

tau_H_inner = cumtrapz(Z_ref(i_Z_150:end),g_ref(i_Z_150:end)*AirComp.BL{'H','M_i'}./(c.R_star*T_ref(i_Z_150:end)));
tau_H_inner = tau_H_inner - tau_H_inner(i_Z_11 - i_Z_150 + 1);

tau_H_outer = cumtrapz(Z_ref(i_Z_150:end),1000*c.phi./D_H.*(T_ref(i_Z_150:end)/Atmos.Z{'11','T_b'}).^(1 + DiffCoeffs{'H','alpha_i'}).*exp(tau_H_inner));
tau_H_outer = tau_H_outer - tau_H_outer(i_Z_11 - i_Z_150 + 1);
tau_H_outer(i_Z_11+1:end) = 0;

n_H_ref = zeros(length(Z_ref),1);
n_H_ref(i_Z_150:end) = (c.n_H_11 - tau_H_outer).*(Atmos.Z{'11','T_b'}./T_ref(i_Z_150:end)).^(1 + DiffCoeffs{'H','alpha_i'}).*exp(-tau_H_inner);
n_H = interp1(Z_ref,n_H_ref,Z);

%% Number Density Interpolation
% Interpolate integral terms for requested altitudes
tau_N2 = interp1(Z_ref,tau_N2_ref,Z);
tau_O = interp1(Z_ref,tau_O_ref,Z);
tau_O2 = interp1(Z_ref,tau_O2_ref,Z);
tau_Ar = interp1(Z_ref,tau_Ar_ref,Z);
tau_He = interp1(Z_ref,tau_He_ref,Z);

% Compute number densities using integral terms
n_N2 = AirComp.BL{'N2','n_i'}*Atmos.Z{'7','T_b'}./T.*exp(-tau_N2);
n_O = AirComp.BL{'O','n_i'}*Atmos.Z{'7','T_b'}./T.*exp(-tau_O);
n_O2 = AirComp.BL{'O2','n_i'}*Atmos.Z{'7','T_b'}./T.*exp(-tau_O2);
n_Ar = AirComp.BL{'Ar','n_i'}*Atmos.Z{'7','T_b'}./T.*exp(-tau_Ar);
n_He = AirComp.BL{'He','n_i'}*Atmos.Z{'7','T_b'}./T.*exp(-tau_He);

%% Helper Functions
% Define height-dependent functions
    function [T] = T_Z(Z)
        L = discretize(Z,Atmos.Z{:,'Z_b'});
        T = zeros(length(Z),1);
        T(L == 1) = Atmos.Z{'7','T_b'} + Atmos.Z{'7','L_Kb'}*(Z(L == 1) - Atmos.Z{'7','Z_b'});
        T(L == 2) = c.T_c + c.A*(1 - ((Z(L == 2) - Atmos.Z{'8','Z_b'})/c.a).^2).^(1/2);
        T(L == 3) = Atmos.Z{'9','T_b'} + Atmos.Z{'9','L_Kb'}*(Z(L == 3) - Atmos.Z{'9','Z_b'});
        T(L == 4 | L == 5) = c.T_inf - (c.T_inf - Atmos.Z{'10','T_b'})*exp(-lambda*zeta_ref(L == 4 | L == 5));
    end

    function [dT_dZ] = dT_dZ_Z(Z)
        L = discretize(Z,Atmos.Z{:,'Z_b'});
        dT_dZ = zeros(length(Z),1);
        dT_dZ(L == 1) = Atmos.Z{'7','L_Kb'};
        dT_dZ(L == 2) = -c.A/c.a*((Z(L == 2) - Atmos.Z{'8','Z_b'})/c.a).*(1 - ((Z(L == 2) - Atmos.Z{'8','Z_b'})/c.a).^2).^(-1/2);
        dT_dZ(L == 3) = Atmos.Z{'9','L_Kb'};
        dT_dZ(L == 4 | L == 5) = lambda.*(c.T_inf - Atmos.Z{'10','T_b'}).*((c.r_0 + Atmos.Z{'10','Z_b'})./(c.r_0 + Z(L == 4 | L == 5))).^2.*exp(-lambda.*zeta_ref(L == 4 | L == 5));
    end

    function [K] = K_Z(Z)
        L = discretize(Z,[86,95,115,1000]);
        K = zeros(length(Z),1);
        K(L == 1) = c.K_7;
        K(L == 2) = c.K_7.*exp(1 - 400./(400 - (Z(L == 2) - 95).^2));
        K(L == 3) = c.K_10;
    end

    function [D] = D_i(G,sum_n)
        D = DiffCoeffs{G,'a_i'}./sum_n.*(T_ref/273.15).^DiffCoeffs{G,'b_i'};
    end

    function [f] = f_i(G,D,M)
        f = g_ref./(c.R_star*T_ref).*(D./(D + K_ref)).*(AirComp.BL{G,'M_i'} + M.*K_ref./D + DiffCoeffs{G,'alpha_i'}*c.R_star./g_ref.*dT_dZ_ref);
    end

    function [v] = v_i(G,Z)
        v = FluxCoeffs{G,'Q_i'}*(Z - FluxCoeffs{G,'U_i'}).^2.*exp(-FluxCoeffs{G,'W_i'}*(Z - FluxCoeffs{G,'U_i'}).^3);
        if strcmp(G,'O')
            v(Z <= 97) = v(Z <= 97) + FluxCoeffs{G,'q_i'}*(FluxCoeffs{G,'u_i'} - Z(Z <= 97)).^2.*exp(-FluxCoeffs{G,'w_i'}*(FluxCoeffs{G,'u_i'} - Z(Z <= 97)).^3);
        end
    end

end
