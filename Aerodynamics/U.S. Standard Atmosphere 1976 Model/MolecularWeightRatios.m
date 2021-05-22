function [MolWeightRatios] = MolecularWeightRatios(u)
%MOLECULARWEIGHTRATIOS Defines molecular weight ratios at geopotential and geometric altitudes.
%
%   Notes:      - In addition to Table 8 data, a Z = 0 and H = 0 data entry was added to facilitate interpolation.
%               - The output structure includes a table for geopotential altitudes and another for geometric altitudes.
%
%   Inputs:     u                       -   universal parameters structure
%
%   Outputs:    MolWeightRatios         -   molecular weight ratios tables
%
%   Author:     William Gravel
%   Created:    03/21/2021
%   Edited:     05/19/2021
%   Purpose:    COESA's U.S. Standard Atmosphere 1976 Model

%% Setup
% Extract universal constants, equations, and tables
e = u.e;

%% Geopotential Altitude-Specific Table (Table 8)
% Define molecular weight ratios
H_ref = [0, 79:0.5:84.5];
Z_ref = e.Z(H_ref/1000)*1000;
M_M_0 = [1.000000, 1.000000, 0.999996, 0.999988, 0.999969, 0.999938, 0.999904, 0.999864, 0.999822, 0.999778, 0.999731, 0.999681, 0.999679];

% Setup table
MolWeightRatios.H = table(H_ref',Z_ref',M_M_0','VariableNames',{'H','Z','M/M_0'});
MolWeightRatios.H.Properties.VariableUnits = ["m'","m",""];
MolWeightRatios.H.Properties.VariableDescriptions = {'Geopotential height','Geometric height','Molecular-weight ratio'};

%% Geometric Altitude-Specific Table (Table 8)
% Define molecular weight ratios
Z_ref = [0, 80:0.5:86];
H_ref = e.H(Z_ref/1000)*1000;
M_M_0 = [1.000000, 1.000000, 0.999996, 0.999989, 0.999971, 0.999941, 0.999909, 0.999870, 0.999829, 0.999786, 0.999741, 0.999694, 0.999641, 0.999579];

% Setup table
MolWeightRatios.Z = table(Z_ref',H_ref',M_M_0','VariableNames',{'Z','H','M/M_0'});
MolWeightRatios.Z.Properties.VariableUnits = ["m","m'",""];
MolWeightRatios.Z.Properties.VariableDescriptions = {'Geometric height','Geopotential height','Molecular-weight ratio'};

end
