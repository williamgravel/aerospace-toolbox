function [DiffCoeffs,FluxCoeffs] = GasSpecies()
%GASSPECIES Define diffusion and flux coefficients for gas species above boundary layer.
%
%   Notes:      - The q_i coefficient is zero for atomic oxygen above Z = 97 km.
%               - The q_i coefficients for all other gases is always zero, and thus u_i and w_i are not defined.
%
%   Inputs:     None
%
%   Outputs:    DiffCoeffs              -   diffusion coefficients table
%               FluxCoeffs              -   flux coefficients table
%
%   Author:     William Gravel
%   Created:    03/21/2021
%   Edited:     05/19/2021
%   Purpose:    COESA's U.S. Standard Atmosphere 1976 Model

%% Diffusion Coefficients (Table 6)
% Define coefficients
G_i = {'N2','O','O2','Ar','He','H'};
alpha_i = [0.00; 0.00; 0.00; 0.00; -0.40; -0.25];
a_i = [NaN; 6.986e20; 4.863e20; 4.487e20; 1.700e21; 3.305e21];
b_i = [NaN; 0.750; 0.750; 0.870; 0.691; 0.500];

% Setup table
DiffCoeffs = table(alpha_i,a_i,b_i,'RowNames',G_i);
DiffCoeffs.Properties.VariableUnits = {'','m^-1*s^-1',''};

%% Flux Coefficients (Table 7)
% Define coefficients
G_i = {'O','O2','Ar','He'};
Q_i = [-5.809644e-4; 1.366212e-4; 9.434079e-5; -2.457369e-4];
q_i = [-3.416248e-3; 0; 0; 0];
U_i = [56.90311; 86.000; 86.000; 86.000];
u_i = [97.0; NaN; NaN; NaN];
W_i = [2.706240e-5; 8.333333e-5; 8.333333e-5; 6.666667e-4];
w_i = [5.008765e-4; NaN; NaN; NaN];

% Setup table
FluxCoeffs = table(Q_i,q_i,U_i,u_i,W_i,w_i,'RowNames',G_i);
FluxCoeffs.Properties.VariableUnits = {'km^-3','km^-3','km','km','km^-3','km^-3'};

end
