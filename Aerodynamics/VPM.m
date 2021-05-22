function [c_l,C_p] = VPM(x_b,y_b,alpha,opts)
%VPM Solves lifting flow equations using Vortex Panel Method.
%   Calculates sectional lift and pressure coefficients given airfoil geometry and flow properties.
%
%   Notes:      - Airfoil coordinates must start at trailing edge and proceed in clockwise fashion.
%
%   Inputs:     x_b                     -   airfoil boundary points x-coordinates
%               y_b                     -   airfoil boundary points y-coordinates
%               alpha                   -   angle of attack
%               opts.AngleUnits         -   requested angle units
%
%   Outputs:    c_l                     -   sectional lift coefficient
%               C_p                     -   pressure coefficient
%
%   Author:     William Gravel
%   Created:    03/23/2021
%   Edited:     04/28/2021

% Define input arguments
arguments
    x_b (:,1) {mustBeReal}
    y_b (:,1) {mustBeReal}
    alpha (1,1) double {mustBeReal}
    opts.AngleUnits (1,1) string {mustBeMember(opts.AngleUnits,{'deg','rad'})} = 'rad'
end

% Validate airfoil coordinates
if any(size(x_b) ~= size(y_b))
    error('Airfoil x- and y-coordinates must be of equal length.')
elseif x_b(1) ~= x_b(end)
    error('Airfoil coordinates must create closed shape.')
end

% Convert angle inputs to radians if in degrees
if strcmp(opts.AngleUnits,'deg')
    alpha = deg2rad(alpha);
end

% Extract chord length from airfoil coordinates
c = max(x_b) - min(x_b);

% Define number of panels from coordinates resolution
N = length(x_b) - 1;

% Compute airfoil panels and first system of equations solutions
x_c = 0.5*(x_b(1:end-1) + x_b(2:end));
y_c = 0.5*(y_b(1:end-1) + y_b(2:end));
S = sqrt((x_b(2:end) - x_b(1:end-1)).^2 + (y_b(2:end) - y_b(1:end-1)).^2);
theta = atan2(y_b(2:end) - y_b(1:end-1),x_b(2:end) - x_b(1:end-1));
RHS = [sin(theta - alpha); 0];

% Define equations for normal and tangent velocity coefficients
A = @(i,j) -(x_c(i) - x_b(j)).*cos(theta(j)) - (y_c(i) - y_b(j)).*sin(theta(j));
B = @(i,j) (x_c(i) - x_b(j)).^2 + (y_c(i) - y_b(j)).^2;
C = @(i,j) sin(theta(i) - theta(j));
D = @(i,j) cos(theta(i) - theta(j));
E = @(i,j) (x_c(i) - x_b(j)).*sin(theta(j)) - (y_c(i) - y_b(j)).*cos(theta(j));
F = @(i,j) log(1 + (S(j).^2 + 2*A(i,j).*S(j))./B(i,j));
G = @(i,j) atan2(E(i,j).*S(j),B(i,j) + A(i,j).*S(j));
P = @(i,j) (x_c(i) - x_b(j)).*sin(theta(i) - 2*theta(j)) + (y_c(i) - y_b(j)).*cos(theta(i) - 2*theta(j));
Q = @(i,j) (x_c(i) - x_b(j)).*cos(theta(i) - 2*theta(j)) - (y_c(i) - y_b(j)).*sin(theta(i) - 2*theta(j));

% Create subscript meshgrid for matrix operations
[i_mesh,j_mesh] = ndgrid(1:N,2:N);

% Compute normal and tangent velocity coefficient matrices
A_n = zeros(N+1);
A_t = zeros(N,N+1);

A_n(1:N,1) = C_n1(1:N,1);
A_n(1:N,N+1) = C_n2(1:N,N);
A_t(1:N,1) = C_t1(1:N,1);
A_t(1:N,N+1) = C_t2(1:N,N);

A_n(1:N,2:N) = C_n1(i_mesh,j_mesh) + C_n2(i_mesh,j_mesh-1);
A_t(1:N,2:N) = C_t1(i_mesh,j_mesh) + C_t2(i_mesh,j_mesh-1);

A_n(N+1,1) = 1;
A_n(N+1,N+1) = 1;

% Compute dimensionless circulation densities along each panel
gamma_prime = A_n\RHS;

% Compute local dimensionless velocities along each panel
V = cos(theta - alpha) + sum(A_t.*gamma_prime',2);

% Calculate total circulation
Gamma = sum(prod([V,S],2));

% Compute coefficient of pressure and coefficient of lift
C_p = 1 - V.^2;
c_l = 2*Gamma/c;

% Define normal and tangent velocity coefficient functions (necessary for i = j checks)
function [C_xx] = C_n2(i,j)
    [C_xx,i,j] = DimCompat(i,j);
    
    i_neq = i(i ~= j);
    j_neq = j(i ~= j);

    C_xx(i == j) = 1;
    C_xx(i ~= j) = D(i_neq,j_neq) + 0.5*Q(i_neq,j_neq).*F(i_neq,j_neq)./S(j_neq) - (A(i_neq,j_neq).*C(i_neq,j_neq) + D(i_neq,j_neq).*E(i_neq,j_neq)).*G(i_neq,j_neq)./S(j_neq);
end

function [C_xx] = C_n1(i,j)
    [C_xx,i,j] = DimCompat(i,j);
    
    i_neq = i(i ~= j);
    j_neq = j(i ~= j);
    
    C_xx(i == j) = -1;
    C_xx(i ~= j) = 0.5*D(i_neq,j_neq).*F(i_neq,j_neq) + C(i_neq,j_neq).*G(i_neq,j_neq) - C_n2(i_neq,j_neq);
end

function [C_xx] = C_t2(i,j)
    [C_xx,i,j] = DimCompat(i,j);
    
    i_neq = i(i ~= j);
    j_neq = j(i ~= j);

    C_xx(i == j) = pi/2;
    C_xx(i ~= j) = C(i_neq,j_neq) + 0.5*P(i_neq,j_neq).*F(i_neq,j_neq)./S(j_neq) + (A(i_neq,j_neq).*D(i_neq,j_neq) - C(i_neq,j_neq).*E(i_neq,j_neq)).*G(i_neq,j_neq)./S(j_neq);
end

function [C_xx] = C_t1(i,j)
    [C_xx,i,j] = DimCompat(i,j);
    
    i_neq = i(i ~= j);
    j_neq = j(i ~= j);

    C_xx(i == j) = pi/2;
    C_xx(i ~= j) = 0.5*C(i_neq,j_neq).*F(i_neq,j_neq) - D(i_neq,j_neq).*G(i_neq,j_neq) - C_t2(i_neq,j_neq);
end

function [C_xx,i,j] = DimCompat(i,j)
    if ~all(size(i) == size(j))
        if length(i) < length(j)
            i = repmat(i,size(j));
        else
            j = repmat(j,size(i));
        end
    end
    if size(i,1) == 1 && size(i,2) ~= 1
        C_xx = zeros(flip(size(i)));
    else
        C_xx = zeros(size(i));
    end
end

end
