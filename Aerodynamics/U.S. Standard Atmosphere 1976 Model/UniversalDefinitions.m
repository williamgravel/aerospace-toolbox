function [u] = UniversalDefinitions()
%UNIVERSALDEFINITIONS Defines adopted constants, fundamental equations, and data tables for model.
%
%   Notes:      - Some constants are initially missing since they need to be calculated using outside functions.
%
%   Inputs:     None
%
%   Outputs:    u                       -   universal parameters structure
%
%   Author:     William Gravel
%   Created:    03/21/2021
%   Edited:     05/21/2021
%   Purpose:    COESA's U.S. Standard Atmosphere 1976 Model

%% Category I Constants (Table 2)
k = 1.380622e-23; % Boltzmann constant [N*m/K]
N_A = 6.022169e26; % Avogadro constant [kmol^-1]
R_star = 8.31432e3; % gas constant [N*m/(kmol*K)]

%% Category II Constants (Table 2)
g_0 = 9.80665; % sea-level acceleration of gravity [m/s^2]
g_0_prime = 9.80665; % standard geopotential meter to geometric height constant [m^2/(s^2*m')]
P_0 = 1.013250e5; % standard sea-level atmospheric pressure [Pa]
r_0 = 6356.766; % effective Earth radius [km]
T_0 = 288.15; % standard sea-level temperature [K]
S = 110.4; % Sutherland constant [K]
beta = 1.458e-6; % constant in expression of dynamic viscosity [kg/(s*m*K^(1/2)]
gamma = 1.400; % ratio of specific heat of air at constant pressure to constant volume
sigma = 3.65e-10; % mean effective collision diameter [m]

%% Category III Constants (Table 2)
T_inf = 1000; % exospheric temperature [K]
phi = 7.2e11; % vertical flux of atomic hydrogen [m^-2*s^-1]

% Temperature definition for elliptical segment (layer 8-9)
T_c = 263.1905; % [K]
A = -76.3232; % [K]
a = -19.9429; % [km]

% Species-specific number densities [atoms/m^3]
n_O_7 = 8.6e16; % for atomic oxygen
n_H_11 = 8.0e10; % for atomic hydrogen

% Eddy-diffusion coefficients [m^2/s]
K_7 = 1.2e2; % for layer 7
K_10 = 0.0; % for layer 10

%% Fundamental Equations
Gamma = g_0/g_0_prime; % geometric to geopotential factor [m'/m]
g_fun = @(Z) g_0*(r_0./(r_0 + Z)).^2; % gravitational accel. as function of geometric altitude [m/s^2]
H_fun = @(Z) Gamma*(r_0*Z./(r_0 + Z)); % geometric to geopotential altitude conversion [km']
Z_fun = @(H) r_0*H./(Gamma*r_0 - H); % geopotential to geometric altitude conversion [km]

%% Unit Conversion Factors (Table 11)
ft_to_m = 3.048e-1;
m_to_ft = 1/ft_to_m;

lb_to_kg = 0.45359237;
kg_to_lb = 1/lb_to_kg;

slug_to_lb = g_0*m_to_ft;
lb_to_slug = 1/slug_to_lb;

K_to_degR = 9/5;
degR_to_K = 1/K_to_degR;

lbf_to_N = lb_to_kg*g_0;
N_to_lbf = 1/lbf_to_N;

%% Parameters Structures
% Create constants structure
c = struct('k',k,'N_A',N_A,'R_star',R_star,'g_0',g_0,'g_0_prime',g_0_prime,'P_0',P_0,'r_0',r_0,...
    'T_0',T_0,'S',S,'beta',beta,'gamma',gamma,'sigma',sigma,'T_inf',T_inf,'T_c',T_c,'A',A,'a',a,...
    'n_O_7',n_O_7,'n_H_11',n_H_11,'K_7',K_7,'K_10',K_10,'phi',phi);

% Create equations structure
e = struct('g',g_fun,'H',H_fun,'Z',Z_fun);

% Create conversion factors structure
f = struct('ft_to_m',ft_to_m,'m_to_ft',m_to_ft,'lb_to_kg',lb_to_kg,'kg_to_lb',kg_to_lb,...
    'slug_to_lb',slug_to_lb,'lb_to_slug',lb_to_slug,'K_to_degR',K_to_degR,'degR_to_K',degR_to_K,...
    'lbf_to_N',lbf_to_N,'N_to_lbf',N_to_lbf);

% Combine structures
u = struct('c',c,'e',e,'f',f);

end
