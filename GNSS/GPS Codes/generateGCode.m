function [G,G_i] = generateGCode(len,stages,taps,phases,opts)
%GENERATEGCODE Generate maximal length sequences.
%   This function is capable of generating maximal length sequences of a given length for a
%   specific number of stages and specific tap locations. This function is also capable of
%   producing a phased version of the m-code if requested.
%
%   Inputs:
%
%      len                     - code length (1x1 double)
%      stages                  - number of shift register stages (1x1 double)
%      taps                    - locations of register taps (1xM double)
%      phases                  - locations of phase selectors (1xN double)
%      opts.CodeType           - code type (1x1 string) {binary/signed}
%
%   Outputs:
%
%      G                       - m-sequence (1xS double)
%      G_i                     - phased m-sequence (1xS double)
%
%   Notes:
%    - Locations of register taps and phase selectors must not exceed number of stages.
%    * Example: A 10-stage shift register may be tapped or phased at any stages between 1 and 10.
%
%    - Codes can be represented either as (a) binary, which includes bits 0 and 1, or (b) signed,
%      which includes positive or negative 1's depending on original binary representation (where
%      0 = unshifted = +1 and 1 = shifted = -1).
%
%   Author:     William Gravel
%   Created:    09/12/2022
%   Edited:     11/26/2022
%
%   See also GENERATECACODE, COMPUTECYCLICCORR, SELECTSVCODE.

arguments
    len (1,1) double {mustBePositive}
    stages (1,1) double {mustBePositive}
    taps (1,:) double {mustBePositive,mustBeLessThanOrEqual(taps,stages)}
    phases (1,:) double {mustBePositive,mustBeLessThanOrEqual(phases,stages)} = [1,1]
    opts.CodeType (1,1) string {mustBeMember(opts.CodeType,["binary","signed"])} = "binary"
end

% Initialize register
reg = ones(1,stages);

% Initialize regular and phased codes
G = zeros(1,len);
G_i = zeros(1,len);

% Iterate over code length
for i = 1:len
    G(i) = reg(end);
    G_i(i) = mod(sum(reg(phases)),2) == 1;
    reg = [mod(sum(reg(taps)),2) == 1, reg(1:end-1)];
end

% Convert binary codes into signed codes
if strcmp(opts.CodeType,"signed")
    G = -2*G + 1;
    G_i = -2*G_i + 1;
end

end
