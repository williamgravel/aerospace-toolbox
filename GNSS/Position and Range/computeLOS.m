function [LOS_enu] = computeLOS(user_ecef,sat_ecef)
arguments
    user_ecef (:,3) double
    sat_ecef (:,3) double
end

LOS_ecef = (sat_ecef - user_ecef)./vecnorm(sat_ecef - user_ecef,2,2);
LOS_llh = ecef2llh(user_ecef);

phi = LOS_llh(:,1);
lambda = LOS_llh(:,2);

e = dot([-sind(lambda), cosd(lambda), zeros(size(phi))],LOS_ecef,2);
n = dot([-sind(phi).*cosd(lambda), -sind(phi).*sind(lambda), cosd(phi)],LOS_ecef,2);
u = dot([cosd(phi).*cosd(lambda), cosd(phi).*sind(lambda), sind(phi)],LOS_ecef,2);

LOS_enu = [e, n, u];

end
