function [az,el,R] = computeAzElRange(user_ecef,sat_ecef)
arguments
    user_ecef (:,3) double
    sat_ecef (:,3) double
end

if size(user_ecef,1) == 1 && size(sat_ecef,1) > 1
    user_ecef = repmat(user_ecef,[size(sat_ecef,1),1]);
elseif size(sat_ecef,1) == 1 && size(user_ecef,1) > 1
    sat_ecef = repmat(sat_ecef,[size(user_ecef,1),1]);
end

LOS_enu = computeLOS(user_ecef,sat_ecef);

az = atan2d(LOS_enu(:,1),LOS_enu(:,2));
el = asind(LOS_enu(:,3)./vecnorm(LOS_enu,2,2));
R = vecnorm(sat_ecef - user_ecef,2,2);

end
