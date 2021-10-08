function [P,rho,N,n_i] = CalcPressure(Z,T,u,t)
%CALCPRESSURE Calculates pressure, density, and total number density for given geometric altitude.
%
%   Inputs:     Z                       -   geometric altitude
%               T                       -   kinetic temperature
%               u                       -   universal parameters structure
%               t                       -   tables structure
%
%   Outputs:    P                       -   pressure
%               rho                     -   density
%               N                       -   total number density
%               n_i                     -   individual number densities
%
%   Author:     William Gravel
%   Created:    03/22/2021
%   Edited:     05/19/2021
%   Purpose:    COESA's U.S. Standard Atmosphere 1976 Model

%% Setup
% Extract universal constants, equations, and tables
c = u.c;
e = u.e;
Atmos = t.Atmos;
AirComp = t.AirComp;
MolWeightRatios = t.MolWeightRatios;

% Initialize vectors
P = zeros(length(Z),1);
rho = zeros(length(Z),1);
N = zeros(length(Z),1);

%% Lower Atmosphere Computations (0 km to 86 km)
% Determine altitude indices at lower atmosphere
is_lower = Z < 86;

if any(is_lower)
    % Convert geometric to geopotential altitude
    H = e.H(Z(is_lower));

    % Discretize height into layer bins
    bins = discretize(H,Atmos.H{:,'H_b'});

    % Determine altitudes with zero temperature gradient
    no_grad = Atmos.H{bins,'L_Mb'} == 0;
    yes_grad = Atmos.H{bins,'L_Mb'} ~= 0;

    % Compute pressure
    P(yes_grad) = Atmos.H{bins(yes_grad),'P_b'}.*(Atmos.H{bins(yes_grad),'T_Mb'}./(Atmos.H{bins(yes_grad),'T_Mb'} + Atmos.H{bins(yes_grad),'L_Mb'}.*(H(yes_grad) - Atmos.H{bins(yes_grad),'H_b'}))).^(c.g_0_prime*1000*c.M_0./(c.R_star*Atmos.H{bins(yes_grad),'L_Mb'}));
    P(no_grad) = Atmos.H{bins(no_grad),'P_b'}.*exp(-c.g_0_prime*1000*c.M_0*(H(no_grad) - Atmos.H{bins(no_grad),'H_b'})./(u.c.R_star*Atmos.H{bins(no_grad),'T_Mb'}));
    
    M = interp1(MolWeightRatios.Z{:,'Z'},MolWeightRatios.Z{:,'M/M_0'},Z(is_lower))*c.M_0;
    
    rho(is_lower) = P(is_lower).*M./(c.R_star*T(is_lower));
    
    N(is_lower) = c.N_A*P(is_lower)./(c.R_star*T(is_lower));
end

%% Upper Atmosphere Computations (86 km to 1000 km)
% Determine altitude indices at upper atmosphere
is_upper = Z >= 86;

if any(is_upper)
    % Calculate number density of individual gas species
    [n_N2,n_O,n_O2,n_Ar,n_He,n_H] = CalcNumberDensity(Z(is_upper),T(is_upper),u,t);

    % Compute total number density
    N(is_upper) = sum([n_N2,n_O,n_O2,n_Ar,n_He,n_H],2);

    % Compute density
    rho(is_upper) = sum([n_N2,n_O,n_O2,n_Ar,n_He,n_H].*AirComp.BL{{'N2','O','O2','Ar','He','H'},'M_i'}',2)/c.N_A;

    % Compute pressure
    P(is_upper) = N(is_upper).*c.k.*T(is_upper);

    % Construct collection of number densities of individual gas species
    n_i = struct('n_N2',n_N2,'n_O',n_O,'n_O2',n_O2,'n_Ar',n_Ar,'n_He',n_He,'n_H',n_H);
end

end
