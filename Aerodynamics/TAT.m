function [a_0,alpha_L_0,c_l,c_m_le,c_m_qc] = TAT(shape,alpha,m,p,opts)
%TAT Solves airfoil aerodynamic properties using Thin Airfoil Theory.
%   Calculates lift slope, zero-lift angle of attack, and sectional lift and moment coefficients for given angle of attack and airfoil geometry.
%
%   Inputs:     shape                   -   airfoil symmetrical/cambered shape
%               alpha                   -   angle of attack
%               m                       -   maximum camber (as percentage of chord)
%               p                       -   position of maximum camber (as tenths of chord)
%               opts.AngleUnits         -   requested angle units
%
%   Outputs:    a_0                     -   lift slope
%               alpha_L_0               -   zero-lift angle of attack
%               c_l                     -   sectional lift coefficient
%               c_m_le                  -   sectional moment coefficient at leading edge
%               c_m_qc                  -   sectional moment coefficient at quarter chord
%
%   Author:     William Gravel - University of Colorado at Boulder
%   Created:    03/25/2021
%   Edited:     04/28/2021

% Define input arguments
arguments
    shape (1,1) string {mustBeMember(shape,{'symmetrical','cambered'})}
    alpha (1,1) double {mustBeReal}
    m (1,1) double {mustBeReal,mustBePositive,mustBeInteger,mustBeInRange(m,0,9)}
    p (1,1) double {mustBeReal,mustBePositive,mustBeInteger,mustBeInRange(p,0,9)}
    opts.AngleUnits (1,1) string {mustBeMember(opts.AngleUnits,{'deg','rad'})} = 'rad'
end

% Convert angle inputs to radians if in degrees
if strcmp(opts.AngleUnits,'deg')
    alpha = deg2rad(alpha);
end


% Define thin airfoil lift slope
a_0 = 2*pi;

% Perform shape-dependent computations
if strcmp(shape,'symmetrical')
    % Compute zero-lift angle of attack
    alpha_L_0 = 0;
elseif strcmp(shape,'cambered')
    % Interpret NACA digits
    m = m/100;
    p = p/10;
    
    % Define coordinate conversion relation   
    x = @(theta) 1/2*(1 - cos(theta));
    
    % Define mean camber line derivative functions
    dz_dx_a = @(theta) 2*m/p^2*(p - x(theta));
    dz_dx_b = @(theta) 2*m/(1 - p)^2*(p - x(theta));

    % Compute angle critical bounds for pre- and post-maximum camber point
    theta_0 = 0;
    theta_1 = acos(1 - 2*p);
    theta_2 = pi;
    
    % Compute zero-lift angle of attack
    alpha_L_0 = -1/pi*(integral(@(theta) dz_dx_a(theta).*(cos(theta) - 1),theta_0,theta_1) + integral(@(theta) dz_dx_b(theta).*(cos(theta) - 1),theta_1,theta_2));
    
    % Compute first non-zero pair of Fourier coefficients
    A_1 = 2/pi*(integral(@(theta) dz_dx_a(theta)*cos(theta),theta_0,theta_1) + integral(@(theta) dz_dx_b(theta)*cos(theta),theta_1,theta_2));
    A_2 = 2/pi*(integral(@(theta) dz_dx_a(theta)*cos(2*theta),theta_0,theta_1) + integral(@(theta) dz_dx_b(theta)*cos(2*theta),theta_1,theta_2));
end

% Determine sectional lift coefficient
c_l = a_0*(alpha - alpha_L_0);

% Calculate sectional moment coefficient at leading edge
if strcmp(shape,'symmetrical')
    c_m_le = -c_l/4;
elseif strcmp(shape,'cambered')
    c_m_le = -(c_l/4 + pi/4*(A_1 - A_2));
end

% Calculate sectional moment coefficient at quarter chord
c_m_qc = c_m_le + c_l/4;

% Convert angles to degrees if requested
if strcmp(opts.AngleUnits,'deg')
    a_0 = deg2rad(a_0);
    alpha_L_0 = rad2deg(alpha_L_0);
end

end
