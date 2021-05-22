function [theta,beta,M] = ThetaBetaM(val,var,opts)
%THETABETAM Solves incomplete oblique shock wave propogation properties.
%   Completes deflection angle, shock wave angle, and Mach number relationship.
%
%   Inputs:     props                   -   incomplete oblique shock wave properties
%               opts.AngleUnits         -   requested angle units
%               opts.Constants          -   custom universal constants
%
%   Outputs:    props                   -   complete oblique shock wave properties
%
%   Author:     William Gravel - University of Colorado at Boulder
%   Created:    04/18/2021
%   Edited:     04/23/2021

% Define input arguments
arguments (Repeating)
    val (1,1) double
    var (1,1) string {mustBeMember(var,{'theta','beta','M'})}
end
arguments
    opts.AngleUnits (1,1) string {mustBeMember(opts.AngleUnits,{'deg','rad'})} = 'deg'
    opts.Constants (1,1) struct = UniversalConstants()
end

% Extract universal constants
gamma = opts.Constants.gamma;

% Check input variables
givens = zeros(3,1);
for i = 1:length(var)
    if strcmp('theta',var{i})
        theta = val{i};
        givens(1) = 1;
    elseif strcmp('beta',var{i})
        beta = val{i};
        givens(2) = 1;
    elseif strcmp('M',var{i})
        M = val{i};
        givens(3) = 1;
    end
end

% Convert angle inputs to degrees if in radians
if strcmp(opts.AngleUnits,'rad')
    if givens(1)
        theta = rad2deg(theta);
    end
    if givens(2)
        beta = rad2deg(beta);
    end
end

% Complete Theta-Beta-M relations
if givens(1) && givens(2) && givens(3)
    error('Too many given variables. No explicit unknown left to solve.')
elseif givens(1) && givens(2)
    theta_lim = atand((5*sind(2*beta))/(5*cosd(2*beta) + 7));
    if theta > theta_lim
        error(sprintf(strcat('Incompatible deflection and shock wave angles.\n',...
            'There is no real Mach number solution for given parameters.')))
    end
    eqn = @(M) atand(2*cotd(beta)*(M^2*sind(beta)^2 - 1)/(M^2*(gamma + cosd(2*beta)) + 2)) - theta;
    M = fzero(eqn,[1,100]);
elseif givens(1) && givens(3)
    beta_max = acosd(((16*M^2)/7 - (4*3^(1/2)*(3*M^4 + 4*M^2 + 20)^(1/2))/7 + 20/7)^(1/2)/(2*M));
    theta_max = atand(2*cotd(beta_max)*(M^2*sind(beta_max)^2 - 1)/(M^2*(gamma + cosd(2*beta_max)) + 2));
    if theta > theta_max
        error(sprintf(strcat('Given deflection angle is larger than maximum possible deflection angle.\n',...
            'There is no straight oblique shock wave solution. A curved bow shock wave is instead formed.')))
    end
    eqn = @(beta) atand(2*cotd(beta)*(M^2*sind(beta)^2 - 1)/(M^2*(gamma + cosd(2*beta)) + 2)) - theta;
    beta_weak = fzero(eqn,[asind(1/M),beta_max]);
    beta_strong = fzero(eqn,[beta_max,90]);
    beta = [beta_weak, beta_strong];
elseif givens(2) && givens(3)
    theta = atand(2*cotd(beta)*(M^2*sind(beta)^2 - 1)/(M^2*(gamma + cosd(2*beta)) + 2));
else
    error('Not enough given variables. Requires two known variables to solve for last unknown.')
end

% Convert angle outputs to radians if requested
if strcmp(opts.AngleUnits,'rad')
    theta = deg2rad(theta);
    beta = deg2rad(beta);
end

end
