function [AirComp,u] = AirComposition(u)
%AIRCOMPOSITION Defines sea-level and boundary layer air composition for present gas species.
%
%   Notes:      - The COESA's model documentation seems to omit the normalization of gas species
%                 fractional shares, resulting in a slightly lower sea-level mean molecular weight.
%
%   Inputs:     u                       -   universal parameters structure
%
%   Outputs:    AirComp                 -   air composition tables
%               u                       -   updated universal parameters structure
%
%   Author:     William Gravel
%   Created:    03/21/2021
%   Edited:     05/19/2021
%   Purpose:    COESA's U.S. Standard Atmosphere 1976 Model

%% Setup
% Extract universal constants, equations, and tables
c = u.c;

%% Sea-Level (H = 0 km) Air Composition (Table 3)
% Define sea-level air composition
G_i = {'N2'; 'O2'; 'Ar'; 'CO2'; 'Ne'; 'He'; 'Kr'; 'Xe'; 'CH4'; 'H2'};
F_i = [0.78084; 0.209476; 0.00934; 0.000314; 0.00001818; 0.00000524; 0.00000114; 0.000000087; 0.000002; 0.0000005];
M_i = [28.0134; 31.9988; 39.948; 44.00995; 20.183; 4.0026; 83.80; 131.30; 16.04303; 2.01594];

% Setup table
AirComp.SL = table(F_i,M_i,prod([F_i,M_i],2),'RowNames',G_i);
AirComp.SL.Properties.VariableNames = {'F_i','M_i','F_i*M_i'};
AirComp.SL.Properties.VariableUnits = {'','kg/kmol','kg/kmol'};
AirComp.SL.Properties.VariableDescriptions = {'Fractional volume','Molecular weight','Concentration product'};

% Update constants
M_0 = sum(F_i.*M_i)/sum(F_i); % sea-level mean molecular weight [kg/kmol]
% M_0 = sum(F_i.*M_i);
u.c.M_0 = M_0;

%% Mesopause Boundary Layer (Z = 86 km) Air Composition (Table 26)
% Define boundary layer air composition
species_bl = {'N2','O2','Ar','He','O', 'H'};
rho_bl = 6.957880e-6; % air density [kg/m^3]
M_O_7 = AirComp.SL{'O2','M_i'}/2; % atomic oxygen molecular weight [kg/kmol]
M_H_7 = AirComp.SL{'H2','M_i'}/2; % atomic hydrogen molecular weight [kg/kmol]
N = (sum(AirComp.SL{species_bl(1:4),'F_i'})*(c.N_A*rho_bl - c.n_O_7*M_O_7))/sum(AirComp.SL{species_bl(1:4),'F_i*M_i'}) + c.n_O_7; % total number density [m^-3]
epsilon = (c.N_A*rho_bl - c.n_O_7*M_O_7)/(N*sum(AirComp.SL{species_bl(1:4),'F_i*M_i'})); % fractional volume factor

F_i_bl = [epsilon*AirComp.SL{species_bl(1:4),'F_i'}; c.n_O_7/N; 0];
M_i_bl = [AirComp.SL{species_bl(1:4),'M_i'}; M_O_7; M_H_7];
n_i_bl = F_i_bl*N;

% Setup table
AirComp.BL = table(F_i_bl,M_i_bl,prod([F_i_bl,M_i_bl],2),n_i_bl,'RowNames',species_bl);
AirComp.BL.Properties.VariableNames = ["F_i'","M_i","F_i'*M_i","n_i"];
AirComp.BL.Properties.VariableUnits = {'','kg/kmol','kg/kmol','m^-3'};
AirComp.BL.Properties.VariableDescriptions = {'Fractional volume','Molecular weight','Concentration product','Number density'};

% Update constants
M_7 = epsilon*sum(AirComp.SL{species_bl(1:4),'F_i*M_i'}) + (c.n_O_7*M_O_7)/N; % boundary layer mean molecular weight [kg/kmol]
u.c.M_7 = M_7;

end
    