function [c_l,c_d] = ShockExpansion(M_inf,alpha,epsilon,shape,opts)
%SHOCKEXPANSION Computes lift and drag coefficients for flow around airfoil.
%   Solves shock expansion problem for flat plate or diamond airfoil given Mach number and angle of attack.
%
%   Inputs:     M_inf                   -   freestream Mach number
%               alpha                   -   angle of attack
%               epsilon                 -   diamond airfoil angle
%               shape                   -   airfoil shape
%               opts.AngleUnits         -   requested angle units
%               opts.Constants          -   custom universal constants
%
%   Outputs:    c_l                     -   sectional lift coefficient
%               c_d                     -   sectional drag coefficient
%
%   Author:     William Gravel - University of Colorado at Boulder
%   Created:    04/24/2021
%   Edited:     05/21/2021

% Define input arguments
arguments
    M_inf (1,1) double
    alpha (1,1) double
    epsilon (1,1) double = 0
    shape (1,1) string {mustBeMember(shape,{'flat','diamond'})} = 'flat'
    opts.AngleUnits (1,1) string {mustBeMember(opts.AngleUnits,{'deg','rad'})} = 'deg'
    opts.Constants (1,1) struct = UniversalConstants()
end

% Extract universal constants
gamma = opts.Constants.gamma;

% Determine freestream properties
[~,nu_inf] = PrandtlMeyer(M_inf,'M','Constants',opts.Constants);

% Determine pressure ratio at upper leading edge
theta_1 = abs(alpha - epsilon);
nu_1 = nu_inf + theta_1;
M_1 = PrandtlMeyer(nu_1,'nu','Constants',opts.Constants);
[~,p_ratio_1_inf] = MachRatios(M_inf,M_1,'Constants',opts.Constants);

% Determine pressure ratio at lower leading edge
theta_3 = alpha + epsilon;
[~,beta_3,~] = ThetaBetaM(M_inf,'M',theta_3,'theta','Constants',opts.Constants);
M_n_inf = M_inf*sind(beta_3(1));
[~,p_ratio_3_inf] = ShockwaveRelations(M_n_inf,'Constants',opts.Constants);
M_n_3 = sqrt((1 + ((gamma - 1)/2)*M_n_inf^2)/(gamma*M_n_inf^2 - (gamma - 1)/2));
M_3 = M_n_3/sind(beta_3(1) - theta_3);
[~,nu_3] = PrandtlMeyer(M_3,'M','Constants',opts.Constants);

% Calculate lift and drag coefficients for flat plate
if strcmpi(shape,'flat')
    c_l = 2/(gamma*M_inf^2)*(p_ratio_3_inf - p_ratio_1_inf)*cosd(alpha);
    c_d = 2/(gamma*M_inf^2)*(p_ratio_3_inf - p_ratio_1_inf)*sind(alpha);
    return    
end

% Determine pressure ratio at upper trailing edge
theta_2 = 2*epsilon;
nu_2 = nu_1 + theta_2;
M_2 = PrandtlMeyer(nu_2,'nu','Constants',opts.Constants);
[~,p_ratio_2_1] = MachRatios(M_1,M_2,'Constants',opts.Constants);
p_ratio_2_inf = p_ratio_2_1*p_ratio_1_inf;

% Determine pressure ratio at lower trailing edge
theta_4 = 2*epsilon;
nu_4 = nu_3 + theta_4;
M_4 = PrandtlMeyer(nu_4,'nu','Constants',opts.Constants);
[~,p_ratio_4_3] = MachRatios(M_3,M_4,'Constants',opts.Constants);
p_ratio_4_inf = p_ratio_4_3*p_ratio_3_inf;

% Calculate normal and axial force coefficients
c_n = 1/(gamma*M_inf^2)*(-p_ratio_1_inf - p_ratio_2_inf + p_ratio_3_inf + p_ratio_4_inf);
c_a = 1/(gamma*M_inf^2)*(p_ratio_1_inf - p_ratio_2_inf + p_ratio_3_inf - p_ratio_4_inf)*tand(epsilon);

% Calculate lift and drag coefficients for diamond airfoil
c_l = c_n*cosd(alpha) - c_a*sind(alpha);
c_d = c_n*sind(alpha) + c_a*cosd(alpha);

end
