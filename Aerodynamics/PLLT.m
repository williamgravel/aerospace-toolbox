function [e,C_L,C_D_i,delta] = PLLT(b,c,a_0,alpha,alpha_L_0,opts)
%PLLT Solves fundamental equation of Prandtl Lifting-Line Theory.
%   Calculates span efficiency factor, lift coefficient, and induced drag coefficient given wing
%   parameters and evaluation options.
%
%   Inputs:     b                       -   wingspan
%               c                       -   root and tip chord length
%               a_0                     -   root and tip lift slope angle of attack [1/rad]
%               alpha                   -   root and tip geometric angle of attack [rad]
%               alpha_L_0               -   root and tip zero-lift angle of attack [rad]
%               opts.EvalNum            -   number of evaluation points
%               opts.EvalPoints         -   angular locations of evaluation points [rad]
%               opts.WingPlanform       -   wing planform shape
%               opts.AngleUnits         -   requested angle units
%
%   Outputs:    e                       -   span efficiency factor
%               C_L                     -   lift coefficient
%               C_D_i                   -   induced drag coefficient
%               delta                   -   induced drag factor
%
%   Author:     William Gravel - University of Colorado at Boulder
%   Created:    03/15/2021
%   Edited:     04/28/2021

% Define input arguments
arguments
    b (1,1) double {mustBeReal,mustBePositive}
    c (1,:) double {mustBeReal,mustBePositive}
    a_0 (1,:) double {mustBeReal}
    alpha (1,:) double {mustBeReal}
    alpha_L_0 (1,:) double {mustBeReal}
    opts.EvalNum (1,1) double {mustBeReal,mustBePositive,mustBeInteger} = 200
    opts.EvalPoints (1,:) double {mustBeReal}
    opts.WingPlanform (1,1) string {mustBeMember(opts.WingPlanform,{'rectangular','elliptical','semi-elliptical'})} = 'rectangular'
    opts.AngleUnits (1,1) string {mustBeMember(opts.AngleUnits,{'deg','rad'})} = 'rad'
end

% Convert angle inputs to radians if in degrees
if strcmp(opts.AngleUnits,'deg')
    a_0 = rad2deg(a_0);
    alpha = deg2rad(alpha);
    alpha_L_0 = deg2rad(alpha_L_0);
    if isfield(opts,'EvalPoints')
        opts.EvalPoints = deg2rad(opts.EvalPoints);
    end
end

% Construct spanwise-varying functions
[c,a_0,alpha,alpha_L_0] = SpanwiseAirfoil(b,c,a_0,alpha,alpha_L_0,opts.WingPlanform);

% Define wing planform geometry
S = integral(c,-b/2,b/2); % planform area
AR = b^2/S; % aspect ratio

% Define evaluation points
if isfield(opts,'EvalNum')
    theta = (1:opts.EvalNum)'*pi/(2*opts.EvalNum);
elseif isfield(opts,'EvalPoints')
    theta = opts.EvalPoints;
end
N = length(theta);
y = b/2*cos(theta);
n = 2*(1:N) - 1;

% Construct grid of evaluation points and Fourier series modes
[theta_grid,n_grid] = ndgrid(theta,n);

% Calculate coefficients of each summation term and populate equation matrix
sys = 4*b./(a_0(y).*c(y)).*sin(n_grid.*theta_grid) + n_grid.*sin(n_grid.*theta_grid)./sin(theta_grid);

% Solve system of equations to find A coefficients
alpha_vec = alpha(y) - alpha_L_0(y);
A = sys\alpha_vec;

% Compute span efficiency factor
delta = sum(n(2:end)'.*(A(2:end)/A(1)).^2);
e = (1 + delta)^(-1);

% Compute lift and induced drag coefficients
C_L = A(1)*pi*AR;
C_D_i = C_L^2/(pi*e*AR);

end
