function [C_p] = CompressibilityCorrection(C_p_0,M_inf,opts)
%COMPRESSIBILITYCORRECTION Corrects pressure coefficients for compressibility.
%
%   Inputs:     C_p_0                   -   incompressible pressure coefficient (Mx1 double)
%               M_inf                   -   freestream Mach number (Mx1 double)
%               opts.Rule               -   compressibility rule to use (1x1 string)
%               opts.RemoveCorrection   -   whether to backwards solve for incompressible pressure (1x1 logical)
%               opts.Constants          -   custom universal constants (1x1 struct)
%
%   Outputs:    C_p                     -   compressible pressure coefficient (Mx1 double)
%
%   Author:     William Gravel - University of Colorado at Boulder
%   Created:    04/20/2021
%   Edited:     04/21/2021

% Define input arguments
arguments
    C_p_0 (:,1) double
    M_inf (:,1) double
    opts.Rule (1,1) string {mustBeMember(opts.Rule,{'Prandtl-Glauert','Karman-Tsien','Laitone'})} = 'Prandtl-Glauert'
    opts.RemoveCorrection (1,1) logical = false
    opts.Constants (1,1) struct = UniversalConstants()
end

% Extract universal constants
gamma = opts.Constants.gamma;

if opts.RemoveCorrection
    % Remove compressibility correction from pressure coefficient
    if strcmp(opts.Rule,'Prandtl-Glauert')
        C_p = C_p_0*sqrt(1 - M_inf^2);
    elseif strcmp(opts.Rule,'Karman-Tsien')
        C_p = -(C_p_0*sqrt(1 - M_inf^2))/((C_p_0*M_inf^2)/(2*(sqrt(1 - M_inf^2) + 1)) - 1);
    elseif strcmp(opts.Rule,'Laitone')
        C_p = -(C_p_0*sqrt(1 - M_inf^2))/((C_p_0*M_inf^2*sqrt(1 - M_inf^2)*((gamma/2 - 1/2)*M_inf^2 + 1))/2 - 1);
    end
else
    % Apply compressibility correction to pressure coefficient
    if strcmp(opts.Rule,'Prandtl-Glauert')
        C_p = C_p_0/sqrt(1 - M_inf^2);
    elseif strcmp(opts.Rule,'Karman-Tsien')
        C_p = C_p_0/(sqrt(1 - M_inf^2) + (M_inf^2/(1 + sqrt(1 - M_inf^2)))*C_p_0/2);
    elseif strcmp(opts.Rule,'Laitone')
        C_p = C_p_0/(sqrt(1 - M_inf^2) + (M_inf^2*(1 + ((gamma - 1)/2)*M_inf^2)/2*sqrt(1 - M_inf^2))*C_p_0);
    end
end

end
