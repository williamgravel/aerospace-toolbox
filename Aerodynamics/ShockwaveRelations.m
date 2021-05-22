function [T_ratio,p_ratio,rho_ratio,T_0_ratio,p_0_ratio,Delta_s] = ShockwaveRelations(M,opts)
%SHOCKWAVERELATIONS Computes thermodynamic ratios across shockwave given flow speed.
%   Outputs downstream to upstream ratios of thermodynamic properties.
%
%   Inputs:     M                       -   Mach number
%               opts.Constants          -   custom universal constants
%
%   Outputs:    T_ratio                 -   temperature ratio
%               p_ratio                 -   pressure ratio
%               rho_ratio               -   density ratio
%               T_0_ratio               -   total temperature ratio
%               p_0_ratio               -   total pressure ratio
%               Delta_s                 -   entropy change
%
%   Author:     William Gravel - University of Colorado at Boulder
%   Created:    04/20/2021
%   Edited:     04/22/2021

% Define input arguments
arguments
    M (1,1) double
    opts.Constants (1,1) struct = UniversalConstants()
end

% Extract universal constants
gamma = opts.Constants.gamma;
c_p = opts.Constants.c_p;
R = opts.Constants.R;

% Calculate shockwave ratios
p_ratio = 1 + 2*gamma/(gamma + 1)*(M^2 - 1);
rho_ratio = (gamma + 1)*M^2/(2 + (gamma - 1)*M^2);
T_ratio = p_ratio*rho_ratio;
T_0_ratio = 1;
Delta_s = c_p*log(T_ratio) - R*log(p_ratio);
p_0_ratio = exp(-Delta_s/R);
end
