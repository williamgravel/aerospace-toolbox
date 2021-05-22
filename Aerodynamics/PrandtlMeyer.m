function [M,nu,mu] = PrandtlMeyer(val,var,opts)
%PRANDTLMEYER Solves incomplete Prandtl-Meyer expansion fan properties.
%   Completes Mach number, Prandtl-Meyer angle, and Mach angle relationships.
%
%   Inputs:     val                     -   given variable value
%               var                     -   given variable name
%               opts.AngleUnits         -   requested angle units
%               opts.Constants          -   custom universal constants
%
%   Outputs:    M                       -   Mach number
%               nu                      -   Prandtl-Meyer angle
%               mu                      -   Mach angle
%
%   Author:     William Gravel - University of Colorado at Boulder
%   Created:    04/18/2021
%   Edited:     04/22/2021

% Define input arguments
arguments
    val (1,1) double
    var (1,1) string {mustBeMember(var,{'M','nu','mu'})} = 'M'
    opts.AngleUnits (1,1) string {mustBeMember(opts.AngleUnits,{'deg','rad'})} = 'deg'
    opts.Constants (1,1) struct = UniversalConstants()
end

% Extract universal constants
gamma = opts.Constants.gamma;

% Convert angle input to degrees if in radians
if strcmp(opts.AngleUnits,'rad') && ismember(var,{'nu','mu'})
    var = rad2deg(var);
end

% Solve Prandtl-Meyer function
nu_fun = @(M) sqrt((gamma + 1)/(gamma - 1))*atand(sqrt((gamma - 1)/(gamma + 1)*(M^2 - 1))) - atand(sqrt(M^2 - 1));

if strcmp(var,'M')
    M = val;
    nu = nu_fun(M);
    mu = asind(1/M);
elseif strcmp(var,'nu')
    nu = val;
    nu_solve = @(M) nu_fun(M) - nu;
    M = fzero(nu_solve,[1,100]);
    mu = asind(1/M);
elseif strcmp(var,'mu')
    mu = val;
    M = 1/sind(mu);
    nu = nu_fun(M);
end

% Convert angle outputs to radians if requested
if strcmp(opts.AngleUnits,'rad')
    nu = deg2rad(nu);
    mu = deg2rad(mu);
end

end
