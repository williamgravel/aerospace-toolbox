function [x,y,x_0,y_t,y_c] = NACA_XXXX(m,p,t,c,N,opts)
%NACA_XXXX Computes x- and y-coordinates of NACA four-digit series airfoil.
%
%   Notes:      - Upper and lower surface coordinates will not equal to reference stations
%                 after rotation transformation is applied
%
%   Inputs:     m                       -   maximum camber (as percentage of chord)
%               p                       -   position of maximum camber (as tenths of chord)
%               t                       -   maximum thickness (as percentage of chord)
%               c                       -   chord length
%               N                       -   number of panels per surface
%               opts.SpacingMethod      -   spacing method for reference stations
%               opts.ClosedTrailingEdge -   whether trailing edge should have zero thickness
%
%   Outputs:    x                       -   whole airfoil x-coordinates
%               y                       -   whole airfoil y-coordinates
%               x_0                     -   reference stations x-coordinates
%               y_t                     -   airfoil thickness y-coordinates
%               y_c                     -   airfoil camberline y-coordinates
%
%   Author:     William Gravel - University of Colorado at Boulder
%   Created:    03/23/2021
%   Edited:     04/25/2021

arguments
    m (1,1) double {mustBeReal,mustBePositive,mustBeInteger,mustBeInRange(m,0,9)}
    p (1,1) double {mustBeReal,mustBePositive,mustBeInteger,mustBeInRange(p,0,9)}
    t (1,1) double {mustBeReal,mustBePositive,mustBeInteger,mustBeInRange(t,0,99)}
    c (1,1) double {mustBeReal,mustBePositive} = 1
    N (1,1) double {mustBeReal,mustBePositive,mustBeInteger} = 100
    opts.SpacingMethod (1,1) string {mustBeMember(opts.SpacingMethod,{'uniform','cosine'})} = 'cosine'
    opts.ClosedTrailingEdge (1,1) logical = false
end

% Interpret NACA digits
m = m/100;
p = p/10;
t = t/100;

% Determine number of points
n = N + 1;

% Create reference stations
if strcmp(opts.SpacingMethod,'uniform')
    x_0 = linspace(1,0,n)';
elseif strcmp(opts.SpacingMethod,'cosine')
    beta = linspace(pi,0,n)';
    x_0 = (1 - cos(beta))/2;
end

% Compute airfoil thickness
if opts.ClosedTrailingEdge
    y_t = t/0.20*(0.2969*x_0.^(1/2) - 0.1260*x_0 - 0.3516*x_0.^2 + 0.2843*x_0.^3 - 0.1036*x_0.^4);
else
    y_t = t/0.20*(0.2969*x_0.^(1/2) - 0.1260*x_0 - 0.3516*x_0.^2 + 0.2843*x_0.^3 - 0.1015*x_0.^4);
end

% Define location indices for forward and aft of max camber
forward = x_0 <= p;
aft = x_0 > p;

% Compute airfoil camber
y_c = zeros(length(x_0),1);
dyc_dx = zeros(length(x_0),1);
if all([m,p] ~= 0)
    y_c(forward) = m/p^2*(2*p*x_0(forward) - x_0(forward).^2);
    y_c(aft) = m/(1 - p)^2*((1 - 2*p) + 2*p*x_0(aft) - x_0(aft).^2);
    dyc_dx(forward) = 2*m/p^2*(p - x_0(forward));
    dyc_dx(aft) = 2*m/(1 - p)^2*(p - x_0(aft));
end
theta = atan(dyc_dx);

% Assemble upper and lower surface coordinates
x_U = x_0 - y_t.*sin(theta);
x_L = x_0 + y_t.*sin(theta);

y_U = y_c + y_t.*cos(theta);
y_L = y_c - y_t.*cos(theta);

% Combine into x- and y-coordinate vectors
x = [x_U; flip(x_L(1:end-1))];
y = [y_U; flip(y_L(1:end-1))];

% Dimensionalize coordinates with chord length
x = x*c;
y = y*c;
x_0 = x_0*c;
y_t = y_t*c;
y_c = y_c*c;

end
