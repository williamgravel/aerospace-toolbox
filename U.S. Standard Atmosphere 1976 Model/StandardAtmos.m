function [Results] = StandardAtmos(alt,opts)
%STANDARDATMOS Computes atmospheric properties of air at given altitudes.
%   Uses COESA's U.S. Standard Atmosphere 1976 model to compute temperature, pressure, density,
%   gravitational acceleration, speed of sound, viscosity, thermal conductivity, and number density
%   of air from 0 km to 1000 km altitude.
%
%   Notes:      - Computed thermal conductivity values do not match COESA 1976 report tables
%
%   Inputs:     alt                     -   altitude
%               opts.ReferenceFrame     -   altitude reference type
%               opts.HeightUnit         -   input height unit
%               opts.UnitSystem         -   output unit system
%               opts.IncludeProps       -   output properties whitelist
%               opts.ExcludeProps       -   output properties blacklist
%               opts.OutputFormat       -   outputs results format
%
%   Outputs:    Results                 -   output results
%
%   Author:     William Gravel
%   Created:    03/02/2021
%   Edited:     05/22/2021
%   Purpose:    COESA's U.S. Standard Atmosphere 1976 Model

%% Argument Validation
arguments
    alt (:,1) double
    opts.ReferenceFrame (1,1) string {mustBeMember(opts.ReferenceFrame,{'Geometric','Geopotential'})} = 'Geometric'
    opts.HeightUnit (1,1) string {mustBeMember(opts.HeightUnit,{'m','km','ft','mi'})} = 'm'
    opts.UnitSystem (1,1) string {mustBeMember(opts.UnitSystem,{'SI','English'})}
    opts.IncludeProps (1,:) string {mustBeMember(opts.IncludeProps,{'Z','H','T','T_M','P','rho','g','N','c_s','mu','nu','k'})} = {}
    opts.ExcludeProps (1,:) string {mustBeMember(opts.ExcludeProps,{'Z','H','T','T_M','P','rho','g','N','c_s','mu','nu','k'})} = {}
    opts.OutputFormat (1,1) string {mustBeMember(opts.OutputFormat,{'table','struct'})} = 'table'
end

opts.HeightUnit = char(opts.HeightUnit);
opts.IncludeProps = convertStringsToChars(opts.IncludeProps);
opts.ExcludeProps = convertStringsToChars(opts.ExcludeProps);

if isfield(opts,'UnitSystem')
    opts.UnitSystem = char(opts.UnitSystem);
else
    if any(strcmp(opts.HeightUnit,{'m','km'}))
        opts.UnitSystem = 'SI';
    elseif any(strcmp(opts.HeightUnit,{'ft','mi'}))
        opts.UnitSystem = 'English';
    end
end

%% Constants, Tables, and Equations
% Define universal definitions structure
u = UniversalDefinitions();

% Define air composition
[AirComp,u] = AirComposition(u);

% Define molecular-weight ratios
[MolWeightRatios] = MolecularWeightRatios(u);

% Define reference levels and gradients of linearly segmented temperature-height profile
[Atmos,u] = AtmosphericLayers(u,MolWeightRatios);

% Define gas species coefficients
[DiffCoeffs,FluxCoeffs] = GasSpecies();

% Combine tables into structure
t = struct('AirComp',AirComp,'MolWeightRatios',MolWeightRatios,'Atmos',Atmos,'DiffCoeffs',DiffCoeffs,'FluxCoeffs',FluxCoeffs);

% Extract structures
c = u.c;
e = u.e;
f = u.f;

%% Altitude Parser
% Interpret altitude input
if strcmp(opts.HeightUnit,'m') % [m] -> [km]
    alt = alt/1000;
elseif strcmp(opts.HeightUnit,'ft') % [ft] -> [km]
    alt = alt*3.048e-4;
elseif strcmp(opts.HeightUnit,'mi') % [mi] -> [km]
    alt = alt*5280*3.048e-4;
end

if strcmp(opts.ReferenceFrame,'Geometric')
    Z = alt;
    H = e.H(Z);
elseif strcmp(opts.ReferenceFrame,'Geopotential')
    H = alt;
    Z = e.Z(H);
end

%% Atmospheric Computations
% Compute temperature (kinetic and molecular) at altitude
[T,T_M] = CalcTemp(Z,u,t);

% Compute pressure, density, and number density at altitude
[P,rho,N] = CalcPressure(Z,T,u,t);

% Compute other atmospheric properties at altitude
c_s = (c.gamma*c.R_star*T_M/c.M_0).^(1/2);
g = e.g(Z);
mu = c.beta*T_M.^(3/2)./(T_M + c.S);
nu = mu./rho;
k = 2.64638e-3*T.^(3/2)./(T + 245.4*10.^(-12./T));

%% Unit Conversions
if strcmp(opts.HeightUnit,'m') % [km] -> [m]
    Z = Z*1000;
    H = H*1000;
elseif strcmp(opts.HeightUnit,'ft') % [km] -> [ft]
    Z = Z/3.048e-4;
    H = H/3.048e-4;
elseif strcmp(opts.HeightUnit,'mi') % [km] -> [mi]
    Z = Z/3.048e-4/5280;
    H = H/3.048e-4/5280;
end

if strcmp(opts.UnitSystem,'English')
    T = T*f.K_to_degR; % [K] -> [degR]
    T_M = T_M*f.K_to_degR; % [K] -> [degR]
    P = P*f.N_to_lbf/(f.m_to_ft)^2; % [N/m^2] -> [lb/ft^2]
    rho = rho*f.kg_to_lb*f.lb_to_slug/(f.m_to_ft)^3; % [kg/m^3] -> [slug/ft^3]
    g = g*f.m_to_ft; % [m/s^2] -> [ft/s^2]
    N = N*(f.m_to_ft)^3; % [m^-3] -> [ft^-3]
    c_s = c_s*f.m_to_ft; % [m/s] -> [ft/s]
    mu = mu*f.kg_to_lb*f.lb_to_slug/f.m_to_ft; % [kg/(m*s)] -> [slug/(ft*s)]
    nu = nu*(f.m_to_ft)^2; % [m^2/s] -> [ft^2/s]
    k = k*f.N_to_lbf/f.K_to_degR; % [W/(m*K)] -> [lb/(s*degR)]
end

%% Results Output
% Setup table
Results = table(Z,H,T,T_M,P,rho,g,N,c_s,mu,nu,k);
if strcmp(opts.UnitSystem,'SI')
    Results.Properties.VariableUnits = {opts.HeightUnit,opts.HeightUnit,'K','K','N/m^2','kg/m^3','m/s^2','m^-3','m/s','kg/(m*s)','m^2/s','W/(m*K)'};
elseif strcmp(opts.UnitSystem,'English')
    Results.Properties.VariableUnits = {opts.HeightUnit,opts.HeightUnit,'degR','degR','lb/ft^2','slug/ft^3','ft/s^2','ft^-3','ft/s','slug/(ft*s)','ft^2/s','lb/(s*degR)'};
end
Results.Properties.VariableDescriptions = {'Geometric altitude','Geopotential altitude','Kinetic temperature','Molecular-scale temperature','Pressure','Air density','Gravitational acceleration','Number density','Speed of sound','Dynamic viscosity','Kinematic viscosity','Thermal coefficient'};

if ~isempty(opts.IncludeProps)
    Results = Results(:,opts.IncludeProps);
elseif ~isempty(opts.ExcludeProps)
    Results = Results(:,setdiff({'Z','H','T','T_M','P','rho','g','N','c_s','mu','nu','k'},opts.ExcludeProps,'stable'));
end

if strcmp(opts.OutputFormat,'struct')
    % Setup struct
    Results = table2struct(Results);
end

end
