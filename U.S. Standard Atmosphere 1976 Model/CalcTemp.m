function [T,T_M] = CalcTemp(Z,u,t)
%CALCTEMP Calculates kinetic and molecular-scale temperature for given geometric altitude.
%
%   Notes:      - Molecular-scale temperature is not defined for heights above Z = 86 km
%
%   Inputs:     Z                       -   geometric altitude
%               u                       -   universal parameters structure
%               t                       -   tables structure
%
%   Outputs:    T                       -   kinetic temperature
%               T_M                     -   molecular-scale temperature
%
%   Author:     William Gravel
%   Created:    03/22/2021
%   Edited:     05/21/2021
%   Purpose:    COESA's U.S. Standard Atmosphere 1976 Model

%% Setup
% Extract universal constants, equations, and tables
c = u.c;
e = u.e;
Atmos = t.Atmos;
MolWeightRatios = t.MolWeightRatios;

%% Temperature Calculations
% Calculate geopotential altitude from geometric altitude
H = e.H(Z);

% Define empty vectors for molecular-scale and kinetic temperatures
T_M = nan(length(Z),1);
T = zeros(length(Z),1);

% Determine thermal region bins for input geometric altitudes
R = discretize(Z,[0,86,1000]);

% Compute geopotential-dependent thermal region
if sum(R == 1) > 0
    L_H = discretize(H(R == 1),Atmos.H{:,'H_b'});
    T_M(R == 1) = Atmos.H{L_H,'T_Mb'} + Atmos.H{L_H,'L_Mb'}.*(H(R == 1) - Atmos.H{L_H,'H_b'});
    T(R == 1) = T_M(R == 1).*interp1(MolWeightRatios.Z{:,'Z'},MolWeightRatios.Z{:,'M/M_0'},Z(R == 1));
end

% Compute geometric-dependent thermal region
if sum(R == 2) > 0
    L_Z = discretize(Z(R == 2),Atmos.Z{:,'Z_b'});
    T(R == 2 & L_Z == 1) = Atmos.Z{'7','T_b'} + Atmos.Z{'7','L_Kb'}.*(Z(R == 2 & L_Z == 1) - Atmos.Z{'7','Z_b'});
    T(R == 2 & L_Z == 2) = c.T_c + c.A*(1 - ((Z(R == 2 & L_Z == 2) - Atmos.Z{'8','Z_b'})/c.a).^2).^(1/2);
    T(R == 2 & L_Z == 3) = Atmos.Z{3,'T_b'} + Atmos.Z{3,'L_Kb'}.*(Z(R == 2 & L_Z == 3) - Atmos.Z{3,'Z_b'});
    T(R == 2 & L_Z > 3) = c.T_inf - (c.T_inf - Atmos.Z{4,'T_b'})*exp(-c.lambda*e.xi(Z(R == 2 & L_Z > 3)));
end

end