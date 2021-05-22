function [T_ratio,p_ratio,rho_ratio,A_ratio] = IsentropicRatios(M,opts)
%ISENTROPICRATIOS Computes thermodynamic ratios for isentropic flow given flow speed.
%   Outputs stagnation to local ratios of thermodynamic properties.
%
%   Notes:      - Assuming isentropic flow and calorically perfect gas
%
%   Inputs:     M                       -   Mach number
%               opts.Constants          -   custom universal constants
%
%   Outputs:    T_ratio                 -   total-local temperature ratio
%               p_ratio                 -   total-local pressure ratio
%               rho_ratio               -   total-local density ratio
%               A_ratio                 -   local-throat nozzle area ratio
%
%   Author:     William Gravel - University of Colorado at Boulder
%   Created:    04/08/2021
%   Edited:     04/21/2021

% Define input arguments
arguments
    M (1,1) double
    opts.Constants (1,1) struct = UniversalConstants()
end

% Extract universal constants
gamma = opts.Constants.gamma;

% Calculate isentropic ratios
T_ratio = 1 + (gamma - 1)/2*M^2;
p_ratio = T_ratio^(gamma/(gamma - 1));
rho_ratio = T_ratio^(1/(gamma - 1));
A_ratio = sqrt(1/M^2*(2/(gamma + 1)*T_ratio)^((gamma + 1)/(gamma - 1)));
end
