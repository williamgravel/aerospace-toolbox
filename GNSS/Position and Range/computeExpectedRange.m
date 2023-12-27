function [R,sat_T_ecef_R] = computeExpectedRange(user_R_ecef_R,sat_R_ecef_R,ephem,t_R,PRN,opts)
arguments
    user_R_ecef_R (:,3) double
    sat_R_ecef_R (:,3) double
    ephem
    t_R (:,2) double
    PRN (1,1) double
    opts.ReferenceFrame (1,1) string = 'gps'
    opts.AbsTol (1,1) double = 1e-6
    opts.PrintIterInfo (1,1) logical = false
end

% Handle differently-sized user and satellite position matrices by reshaping smaller matrix
% (note: this only works if either matrix has a single row of data OR both matrix row counts agree)
if size(user_R_ecef_R,1) == 1 && size(sat_R_ecef_R,1) > 1
    user_R_ecef_R = repmat(user_R_ecef_R,[size(sat_R_ecef_R,1),1]);
elseif size(sat_R_ecef_R,1) == 1 && size(user_R_ecef_R,1) > 1
    sat_R_ecef_R = repmat(sat_R_ecef_R,[size(user_R_ecef_R,1),1]);
elseif (size(user_R_ecef_R,1) > 1 && size(sat_R_ecef_R,1) > 1) && size(user_R_ecef_R,1) ~= size(sat_R_ecef_R,1)
    error('MATLAB:sizeDimensionsMustMatch','Arrays have incompatible sizes for this operation')
end

% t_R = repmat(t_R,[size(user_R_ecef_R,1),1]);

% Load parameters for ellipsoid model of choice
model = TRF(opts.ReferenceFrame);

% Initialize iteration count and residual value
i = 1;
res = Inf;

% Compute "a priori" geometric range
R = vecnorm(sat_R_ecef_R - user_R_ecef_R,2,2);

while res > opts.AbsTol
    % Keep previous range values in memory
    R_old = R;

    % Calculate transmission time based on reception time and transit time
    t_T = [t_R(:,1) - R/model.c, t_R(:,2)];

    % Compute satellite position at transmission time
    sat_T_ecef_T = eph2pos(ephem,t_T,PRN);

    % Rotate satellite position from ECEF @ t_T to ECEF @ t_R (due to Earth's rotation rate)
    phi = model.Omega_dot*(t_R(:,1) - t_T(:,1));
    sat_T_ecef_R = [dot([cos(phi),sin(phi),zeros(size(phi))],sat_T_ecef_T,2),...
        dot([-sin(phi),cos(phi),zeros(size(phi))],sat_T_ecef_T,2),...
        dot([zeros(size(phi)),zeros(size(phi)),ones(size(phi))],sat_T_ecef_T,2)];

    % Compute new geometric range
    R = vecnorm(sat_T_ecef_R - user_R_ecef_R,2,2);

    % Calculate range residuals from last iteration to current iteration
    res = max(abs(R - R_old));

    % Print iteration/residual info if requested
    if opts.PrintIterInfo
        fprintf('Iteration #%.f: diff = %.2e [m]\n',i,res)
    end

    % Increment iteration count
    i = i + 1;
end

end
