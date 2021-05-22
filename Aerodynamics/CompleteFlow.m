function [flow] = CompleteFlow(flow)
%COMPLETEFLOW Completes thermodynamic properties of given partial flow.
%
%   Inputs:     flow              -   flow properties
%
%   Outputs:    flow              -   updated flow properties
%
%   Author:     William Gravel - University of Colorado at Boulder
%   Created:    04/09/2021
%   Edited:     05/21/2021

% Define input arguments
arguments
    flow (1,1) struct
end

% Determine missing thermal property
if isfield(flow,'T') && ~isfield(flow,'a')
    flow.a = sqrt(gamma*R*flow.T);
elseif ~isfield(flow,'T') && isfield(flow,'a')
    flow.T = flow.a^2/(gamma*R);
elseif ~isfield(flow,'T') && ~isfield(flow,'a')
    warning('Not enough information. Missing either temperature or speed of sound.')
end

% Determine missing speed property
if isfield(flow,'M') && ~isfield(flow,'u')
    flow.u = flow.M*flow.a;
elseif ~isfield(flow,'M') && isfield(flow,'u')
    flow.M = flow.u/flow.a;
elseif ~isfield(flow,'M') && ~isfield(flow,'u')
    warning('Not enough information. Missing either Mach number or flow velocity.')
end

% Determine missing fluid property
if isfield(flow,'p') && ~isfield(flow,'rho')
    flow.rho = flow.p/(R*flow.T);
elseif ~isfield(flow,'p') && isfield(flow,'rho')
    flow.p = flow.rho*R*flow.T;
elseif ~isfield(flow,'p') && ~isfield(flow,'rho')
    warning('Not enough information. Missing either pressure or density.')
end

end
