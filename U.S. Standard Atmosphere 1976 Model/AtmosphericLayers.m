function [Atmos,u] = AtmosphericLayers(u,MolWeightRatios)
%ATMOSPHERICLAYERS Defines base altitudes, temperature gradients, and initial temperatures and pressures for atmospheric layers.
%
%   Inputs:     u                       -   universal parameters structure
%               MolWeightRatios         -   molecular weight ratios tables
%
%   Outputs:    Atmos                   -   atmospheric layers tables
%               u                       -   updated universal parameters structure
%
%   Author:     William Gravel
%   Created:    03/21/2021
%   Edited:     05/19/2021
%   Purpose:    COESA's U.S. Standard Atmosphere 1976 Model

%% Setup
% Extract universal constants, equations, and tables
c = u.c;
e = u.e;

%% Molecular-Scale Temperatures at Geopotential-Based Atmospheric Layers (Table 4)
% Define layer-specific base altitudes and temperature gradients
H_b = [0; 11; 20; 32; 47; 51; 71; e.H(86)];
L_Mb = [-6.5; 0.0; 1.0; 2.8; 0.0; -2.8; -2.0; NaN];

% Compute layer-specific initial temperatures
T_Mb = [c.T_0; c.T_0 + cumsum(L_Mb(1:end-1).*(H_b(2:end) - H_b(1:end-1)))];

% Compute layer-specific initial pressures
P_b = c.P_0;
for i = 1:7
    if L_Mb(i) ~= 0
        P_b(i+1) = P_b(i)*(T_Mb(i)/(T_Mb(i) + L_Mb(i)*(H_b(i+1) - H_b(i))))^(c.g_0_prime*1000*c.M_0/(c.R_star*L_Mb(i)));
    else
        P_b(i+1) = P_b(i)*exp(-c.g_0_prime*1000*c.M_0*(H_b(i+1) - H_b(i))/(c.R_star*T_Mb(i)));
    end
end
P_b = P_b';

% Setup table
Atmos.H = table(H_b,L_Mb,T_Mb,P_b,'RowNames',string(0:7));
Atmos.H.Properties.VariableUnits = ["km'","K/km'","K","Pa"];
Atmos.H.Properties.VariableDescriptions = {'Geopotential height','Molecular-scale temperature gradient','Molecular-scale temperature','Pressure'};

%% Kinetic Temperatures at Geometric-Based Atmospheric Layers (Table 5)
% Define layer-specific base altitudes and temperature gradients
Z_b = [86; 91; 110; 120; 500; 1000];
L_Kb = [0.0; NaN; 12.0; NaN; NaN; NaN];

% Compute layer-specific initial temperatures
T_b(1) = T_Mb(end)*MolWeightRatios.Z{end,'M/M_0'};
T_b(2) = T_b(1) + L_Kb(1)*(Z_b(2) - Z_b(1));
T_b(3) = c.T_c + c.A*(1 - ((Z_b(3) - Z_b(2))/c.a)^2)^(1/2);
T_b(4) = T_b(3) + L_Kb(3)*(Z_b(4) - Z_b(3));

lambda = L_Kb(3)/(c.T_inf - T_b(4));
xi = @(Z) (Z - Z_b(4))*(c.r_0 + Z_b(4))./(c.r_0 + Z);

T_b(5) = c.T_inf - (c.T_inf - T_b(4))*exp(-lambda*xi(Z_b(5)));
T_b(6) = c.T_inf - (c.T_inf - T_b(4))*exp(-lambda*xi(Z_b(6)));
T_b = T_b';

% Setup table
Atmos.Z = table(Z_b,L_Kb,T_b,'RowNames',string(7:12));
Atmos.Z.Properties.VariableUnits = ["km","K/km","K"];
Atmos.Z.Properties.VariableDescriptions = {'Geometric height','Kinetic-temperature gradient','Kinetic temperature'};

% Update constants and equations
u.c.lambda = lambda;
u.e.xi = xi;

end
