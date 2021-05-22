function [T_ratio,p_ratio,rho_ratio] = MachRatios(M_1,M_2,opts)
%MACHRATIOS Computes thermodynamic ratios for isentropic flow given flow speed change.
%   Outputs downstream to upstream ratios of thermodynamic properties.
%
%   Notes:      - Assuming isentropic flow and calorically perfect gas
%
%   Inputs:     M_1                     -   upstream Mach number
%               M_2                     -   downstream Mach number
%               opts.Constants          -   custom universal constants
%
%   Outputs:    T_ratio                 -   temperature ratio
%               p_ratio                 -   pressure ratio
%               rho_ratio               -   density ratio
%
%   Author:     William Gravel - University of Colorado at Boulder
%   Created:    04/24/2021
%   Edited:     04/24/2021

% Define input arguments
arguments
    M_1 (1,1) double
    M_2 (1,1) double
    opts.Constants (1,1) struct = UniversalConstants()
end

% Extract universal constants
gamma = opts.Constants.gamma;

% Calculate isentropic ratios
T_ratio = (1 + ((gamma - 1)/2)*M_1^2)/(1 + ((gamma - 1)/2)*M_2^2);
p_ratio = T_ratio^(gamma/(gamma - 1));
rho_ratio = T_ratio^(1/(gamma - 1));
end
