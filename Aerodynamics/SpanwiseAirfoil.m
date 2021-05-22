function [c_y,a_0_y,alpha_y,alpha_L_0_y] = SpanwiseAirfoil(b,c,a_0,alpha,alpha_L_0,shape)
%SPANWISEAIRFOIL Constructs spanwise-varying functions for wing and airfoil properties.
%   Outputs functions for chord length, lift slope, geometric angle of attack, and zero-lift angle of attack.
%
%   Notes:      - Spanwise varying quantities may be input as root and tip vector or single constant scalar
%               - Assumes linearly varying quantities from root to tip
%
%   Inputs:     b                       -   wingspan
%               c                       -   root and tip chord length
%               a_0                     -   root and tip lift slope angle of attack
%               alpha                   -   root and tip geometric angle of attack
%               alpha_L_0               -   root and tip zero-lift angle of attack
%               shape                   -   wing planform shape
%
%   Outputs:    c_y                     -   spanwise chord length function
%               a_0_y                   -   spanwise lift slope function
%               alpha_y                 -   spanwise geometric angle of attack function
%               alpha_L_0_y             -   short variable description
%
%   Author:     William Gravel - University of Colorado at Boulder
%   Created:    04/25/2021
%   Edited:     04/27/2021

% Define input arguments
arguments
    b (1,1) double {mustBeReal,mustBePositive}
    c (1,:) double {mustBeReal,mustBePositive}
    a_0 (1,:) double {mustBeReal}
    alpha (1,:) double {mustBeReal}
    alpha_L_0 (1,:) double {mustBeReal}
    shape (1,1) string {mustBeMember(shape,{'rectangular','elliptical','semi-elliptical'})}
end

% Define chord length function
if strcmp(shape,'rectangular')
    if length(c) == 1
        c_y = @(y) c*ones(size(y));
    elseif length(c) == 2
        c_y = @(y) c(1) - (c(2) - c(1))*abs(y)/(b/2);
    else
        error('Invalid chord length input argument size.')
    end
elseif strcmp(shape,'elliptical')
    c_y = @(y) sqrt((1 - y.^2/(b/2)^2)*c(1)^2)*2;
elseif strcmp(shape,'semi-elliptical')
    c_y = @(y) sqrt((1 - y.^2/(b/2)^2)*c(1)^2);
end

% Define lift slope function
if length(a_0) == 1
    a_0_y = @(y) a_0*ones(size(y));
elseif length(a_0) == 2
    a_0_y = @(y) a_0(1) - (a_0(2) - a_0(1))*abs(y)/(b/2);
else
    error('Invalid lift slope input argument size.')
end

% Define geometric angle of attack function
if length(alpha) == 1
    alpha_y = @(y) alpha*ones(size(y));
elseif length(alpha) == 2
    alpha_y = @(y) alpha(1) - (alpha(2) - alpha(1))*abs(y)/(b/2);
else
    error('Invalid geometric angle of attack input argument size.')
end

% Define zero-lift angle of attack function
if length(alpha_L_0) == 1
    alpha_L_0_y = @(y) alpha_L_0*ones(size(y));
elseif length(alpha_L_0) == 2
    alpha_L_0_y = @(y) alpha_L_0(1) - (alpha_L_0(2) - alpha_L_0(1))*abs(y)/(b/2);
else
    error('Invalid zero-lift angle of attack input argument size.')
end

end
